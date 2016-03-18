import sys
import time
from pb_model import Constraint
import collections
from pprint import pprint as pp

class Instance(object):
    def __init__(self, lines, pb_model, max_bp):
        self.pb_model = pb_model
        self.max_bp = max_bp   # Maximum permitted number of blocking pairs

        self.read_lines(lines)
        self.presolve()

        # self.rplace[i][j]==1 <-> resident i gets his j^th choice
        # self.r_unassigned[i]==1 <-> resident i is not assigned to a hospital
        self.rplace = [None] * self.nres
        self.r_unassigned = [None] * self.nres
        for i in [c[0] for c in self.couples] + self.singles:
            self.rplace[i] = [self.pb_model.create_var("res{}-{}".format(i, j))
                                    for j in range(len(self.rpref[i]))]
            self.r_unassigned[i] = self.pb_model.create_var("res{}_unassigned".format(i))
            self.pb_model.add_exactly_one_constr(self.rplace[i] + [self.r_unassigned[i]],
                    "Resident assigned exactly once, or unassigned")

        for i, j in self.couples:
            # Identical vars for second member of each couple
            self.rplace[j] = self.rplace[i]
            self.r_unassigned[j] = self.r_unassigned[i]

        # self.hplace[i][j]==1 <-> hospital i gets its j^th choice
        self.hplace = []
        for i in range(self.nhosp):
            self.hplace.append([self.pb_model.create_var("hosp{}-{}".format(i, j))
                                    for j in range(len(self.hpref[i]))])
            if len(self.hpref[i]) > self.hosp_cap[i]:  # If hospital has more prefs than capacity
                self.pb_model.add_sum_leq_constr(self.hplace[-1], self.hosp_cap[i],
                        "Hospital capacity")

        # Chosen hosp prefs match chosen res prefs
        # TODO: re-use res pref vars in hospital lists for single residents
        for i, prefs in enumerate(self.hpref):
            for pos, res in enumerate(prefs):
                self.pb_model.add_constr(
                        Constraint([(1, self.hplace[i][pos])] +
                                   [(-1, self.rplace[res][k]) for k in self.rrank(res, i)],
                                   "=", 0, "Hosp pref matches res prefs"))

        # Add stability
        self.bp_vars = []   # Blocking pair vars

        for i in self.singles:
            self.add_type1(i)

        for i, j in self.couples:
            self.add_type2(i, j)
            self.add_type2(j, i)
            self.add_type3(i, j)

        self.pb_model.add_sum_leq_constr(self.bp_vars, self.max_bp, "Max permitted number of blocking pairs")

        # Add objective: two points for each matched couple and one point for each matched single
        obj_terms = ([(2, v) for i,_ in self.couples for v in self.rplace[i]] +
                     [(1, v) for i in self.singles for v in self.rplace[i]])
        self.pb_model.add_objective(obj_terms)

    def read_lines(self, lines):
        self.nres = int(lines[0])
        self.nhosp = int(lines[1])
        self.ncoup = int(lines[2])
        self.nsingle = self.nres - 2 * self.ncoup
        self.first_single = 2 * self.ncoup
        self.couples = [(i*2, i*2+1) for i in range(self.ncoup)]
        self.singles = range(self.first_single, self.nres)
        self.npost = int(lines[3])
        lines = lines[9:]

        # Pref list for each resident
        self.rpref = [[int(x) for x in line.split()[1:]] for line in lines[:self.nres]]
        
        # Pref list for each hospital
        self.hpref = [[int(x) for x in line.split()[2:]] for line in lines[self.nres:self.nres+self.nhosp]]

        # Hospital capacities
        self.hosp_cap = [int(line.split()[1]) for line in lines[self.nres:self.nres+self.nhosp]]

    def presolve(self):
        """Reduce the size of the problem by removing some preferences that,
           if chosen, would result in too many blocking pairs.
        """
        q = collections.deque()  # Queue of residents to process
        in_q = [False] * self.nres
        for i in [r1 for r1, r2 in self.couples] + self.singles:
            q.append(i)
            in_q[i] = True

        while q:
            i = q.popleft()
            in_q[i] = False
            # num_bp counts the number of pairs that must be blocking pairs
            # if i doesn't get this or a better choice,
            # due to the low position of i on the hospital's preference list
            num_bp = 0
            if self.is_single(i):
                for j, hosp in enumerate(self.rpref[i][:-1]):
                    if self.hrank(hosp, i) < self.hosp_cap[hosp]:
                        num_bp += 1
                    if num_bp > self.max_bp:
                        #print "Cutting down pref list of res {} after pos {}".format(i, j)
                        for hosp in self.rpref[i][j+1:]:
                            self.hpref[hosp].remove(i)
                            if len(self.hpref[hosp]) >= self.hosp_cap[hosp]:
                                r = self.hpref[hosp][self.hosp_cap[hosp] - 1]
                                if not in_q[r]:
                                    if not self.is_single(r) and r % 2 == 1:
                                        r -= 1  # Always add the first member of a couple
                                    q.append(r)
                                    in_q[r] = True
                        self.rpref[i] = self.rpref[i][:j+1]  # Trim resident's preferences
                        break
            else:
                for j, (hosp1, hosp2) in enumerate(zip(self.rpref[i][:-1], self.rpref[i+1][:-1])):
                    hosp1_cap = self.hosp_cap[hosp1]
                    if hosp1 == hosp2:
                        # Both residents in couples are in hosp1's top-<capacity> choices
                        if self.hrank(hosp1, i) < hosp1_cap and self.hrank(hosp1, i+1) < hosp1_cap:
                            num_bp += 1
                    else:
                        hosp2_cap = self.hosp_cap[hosp2]
                        if self.hrank(hosp1, i) < hosp1_cap and self.hrank(hosp2, i+1) < hosp2_cap:
                            num_bp += 1

                    if num_bp > self.max_bp:
                        # remove residents in this couple from preference lists of hospitals
                        # where they can no longer appear
                        for res in [i, i+1]:
                            # h_remove is the set of hospitals from whose pref lists we can remove resident res
                            h_remove = set(self.rpref[res])
                            for hosp in self.rpref[res][:j+1]: h_remove.discard(hosp)
                            for hosp in h_remove:
                                self.hpref[hosp].remove(res)
                                if len(self.hpref[hosp]) >= self.hosp_cap[hosp]:
                                    r = self.hpref[hosp][self.hosp_cap[hosp] - 1]
                                    if not in_q[r]:
                                        if not self.is_single(r) and r % 2 == 1:
                                            r -= 1  # Always add the first member of a couple
                                        q.append(r)
                                        in_q[r] = True

                            # Trim resident's preferences
                            self.rpref[res] = self.rpref[res][:j+1] 
                        break
#        pp (self.rpref[:self.first_single])
#        print
#        pp (self.rpref[self.first_single:])
#        sys.exit(1)
           
    def write(self, quiet):
        self.pb_model.write(quiet)

    def add_type1(self, i):
        for j, h in enumerate(self.rpref[i]):
            # If hospital h would like to take resident i, hosp_has_space_var must take the value 1.
            # Otherwise, this variable can take any value.
            hosp_has_space_var = self.pb_model.create_var("type1_hosp_space-{}-{}-{}".format(i, h, self.hrank(h, i)))
            self.enforce_hosp_space_var(hosp_has_space_var, h, self.hosp_cap[h], self.hrank(h, i),
                    "Hosp has space var has correct value (type 1)")
            var = self.pb_model.create_var("type1-{}-{}".format(i, j))
            self.bp_vars.append(var)
            self.pb_model.add_constr(Constraint([(1, var), (-1, hosp_has_space_var)] +
                                      [(1, v) for v in self.rplace[i][:j+1]], ">=", 0, "Type 1 stability"))

    def add_type2(self, i, partner):
        for j, h in enumerate(self.rpref[i]):
            hrank_of_res = self.hrank(h, i)
            # If hospital h isn't filled to capacity by i's partner and residents preferred to i,
            # then hosp_has_space_var must take the value 1.
            hosp_has_space_var = self.pb_model.create_var("type2_hosp_space-{}-{}-{}".format(i, h, hrank_of_res))
            hplace_vars = self.hplace[h][:hrank_of_res + 1]
            try:
                hrank_of_partner = self.hrank(h, partner) 
                if hrank_of_partner > hrank_of_res:
                    hplace_vars.append(self.hplace[h][hrank_of_partner])
            except ValueError:
                pass   # Hospital doesn't rank partner
            self.pb_model.add_constr(Constraint([(self.hosp_cap[h], hosp_has_space_var)] + 
                                      [(1, v) for v in hplace_vars],
                                      ">=", self.hosp_cap[h], "Hosp has space var has correct value (type 2)"))
            var = self.pb_model.create_var("type2-{}-{}".format(i, j))
            self.bp_vars.append(var)
            # A list of ranks worse than j for resident i, such that i's partner gets the same hospital as he does at rank j
            worse_ranks_with_same_partner_hosp = [
                    idx for idx in range(j+1, len(self.rpref[i])) if self.rpref[partner][j]==self.rpref[partner][idx]]
            self.pb_model.add_constr(Constraint([(1, var), (-1, hosp_has_space_var)] +
                              [(-1, self.rplace[i][idx]) for idx in worse_ranks_with_same_partner_hosp],
                              ">=", -1, "Type 2 stability"))
        
    def add_type3(self, i, partner):
        for j, (h, h2) in enumerate(zip(self.rpref[i], self.rpref[partner])):
            hrank_of_res = self.hrank(h, i)
            hrank_of_partner = self.hrank(h2, partner) 
            if h != h2:
                self.add_type3a(i, partner, j, h, h2, hrank_of_res, hrank_of_partner)
            elif hrank_of_res < hrank_of_partner:
                self.add_type3bcd(i, partner, j, h, hrank_of_res, hrank_of_partner)
            else:
                self.add_type3bcd(partner, i, j, h, hrank_of_partner, hrank_of_res)

    def add_type3a(self, i, partner, j, h, h2, hrank_of_res, hrank_of_partner):
        hosp1_has_space_var = self.pb_model.create_var("hosp1_space-{}-{}".format(h, hrank_of_res))
        self.enforce_hosp_space_var(hosp1_has_space_var, h, self.hosp_cap[h], hrank_of_res,
                "Hosp 1 has space var has correct value (type 3a)")

        hosp2_has_space_var = self.pb_model.create_var("hosp2_space-{}-{}".format(h2, hrank_of_partner))
        self.enforce_hosp_space_var(hosp2_has_space_var, h2, self.hosp_cap[h2], hrank_of_partner,
                "Hosp 2 has space var has correct value (type 3a)")

        var = self.pb_model.create_var("type3a-{}-{}".format(i, j))
        self.bp_vars.append(var)
        self.pb_model.add_constr(Constraint([(1, var), (-1, hosp1_has_space_var), (-1, hosp2_has_space_var)] +
                      [(-1, self.rplace[i][idx]) for idx in range(j+1, len(self.rpref[i]))
                            if h != self.rpref[i][idx] and h2 != self.rpref[partner][idx]] +
                      [(-1, self.r_unassigned[i])],
                      ">=", -2, "Type 3a stability"))
        
    def add_type3bcd(self, i, partner, j, h, hrank_of_res, hrank_of_partner):
        # The hospital should have no more than capacity-2 spaces used up to resident i,
        # and no more than capacity-1 spaces used up to resident partner
        hosp_has_space_part_1_var = self.pb_model.create_var("hosp(3bcd)_space-{}-{}-part-1".format(h, hrank_of_res))
        self.enforce_hosp_space_var(hosp_has_space_part_1_var, h, self.hosp_cap[h] - 1, hrank_of_res,
                "Hosp (3bcd) has space var part 1 has correct value")

        hosp_has_space_part_2_var = self.pb_model.create_var("hosp(3bcd)_space-{}-{}-part-2".format(h, hrank_of_partner))
        self.enforce_hosp_space_var(hosp_has_space_part_2_var, h, self.hosp_cap[h], hrank_of_partner,
                "Hosp (3bcd) has space var part 2 has correct value")

        # If there's space in both hospitals, then hosp_has_space_var must take value 1.
        hosp_has_space_var = self.pb_model.create_var("hosp(3bcd)_space-{}-{}".format(h, hrank_of_res))
        self.pb_model.add_constr(Constraint([(1, hosp_has_space_var),
                                             (-1, hosp_has_space_part_1_var),
                                             (-1, hosp_has_space_part_2_var)
                                             ], ">=", -1, "Hospital has space (3bcd)"))

        var = self.pb_model.create_var("type3bcd-{}-{}".format(i, j))
        self.bp_vars.append(var)
        self.pb_model.add_constr(Constraint([(1, var), (-1, hosp_has_space_var), (-1, self.r_unassigned[i])] +
                          [(-1, self.rplace[i][idx]) for idx in range(j+1, len(self.rpref[i]))
                                if h != self.rpref[i][idx] and h != self.rpref[partner][idx]],
                          ">=", -1, "Type 3bcd stability"))
        
    def enforce_hosp_space_var(self, v, h, n_res, up_to_pos, constraint_name):
        """Enforce the condition that if hospital h has fewer than n_res residents assigned in positions
           <= up_to_pos in h's preference list, then variable v must take the value 1.
           (Otherwise, v can take any value.)
        """
        self.pb_model.add_constr(Constraint(
                [(n_res, v)] + [(1, v) for v in self.hplace[h][:up_to_pos + 1]], ">=", n_res, constraint_name))
        
    def is_single(self, res):
        return res >= self.first_single

    def rrank(self, r, h):
        """What ranks does resident r give hospital h (as an array)?
           Note that a resident in a couple may rank a hospital more than once
        """
        return [i for i, hosp in enumerate(self.rpref[r]) if hosp == h]
            
    def hrank(self, h, r):
        "What rank does hospital h give resident r?"
        return self.hpref[h].index(r)

