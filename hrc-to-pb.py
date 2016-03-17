import sys
from pprint import pprint as pp

from pb_model import Constraint, PBModel

class Instance(object):
    def __init__(self, lines, pb_model, max_bp):
        self.pb_model = pb_model
        self.max_bp = max_bp   # Maximum permitted number of blocking pairs

        self.read_lines(lines)

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
            # Hospital capacity constraint
            self.pb_model.add_sum_leq_constr(self.hplace[-1], self.hosp_cap[i], "Hospital capacity")

        # Chosen hosp prefs match chosen res prefs
        # TODO: re-use res pref vars in hospital lists for single residents
        for i, prefs in enumerate(self.hpref):
            for pos, res in enumerate(prefs):
                self.pb_model.add_constr(
                        Constraint([(1, self.hplace[i][pos])] +
                                   [(-1, self.rplace[res][k]) for k in self.rrank(res, i)],
                                   "=", 0, "Hosp pref matches res prefs"))

        self.add_stability()
        self.add_objective()

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

    def write(self, quiet):
        self.pb_model.write_model_size_comment()
        if not quiet:
            self.pb_model.show_var_names()
            self.pb_model.show_objective()
        self.pb_model.write_model(quiet)

    def add_objective(self):
        obj_terms = []
        for i in range(self.ncoup):
            for v in self.rplace[i*2]:
                obj_terms.append((2, v))
        for i in range(self.nsingle):
            for v in self.rplace[self.first_single + i]:
                obj_terms.append((1, v))
        self.pb_model.add_objective(obj_terms)
        
    def add_stability(self):
        self.bp_vars = []   # Blocking pair vars

        # TYPE 1
        for i in self.singles:
            for j, h in enumerate(self.rpref[i]):
                hosp_has_space_var = self.pb_model.create_var("hosp_space-{}-{}".format(h, self.hrank(h, i)))
                self.pb_model.add_constr(Constraint([(self.hosp_cap[h], hosp_has_space_var)] + 
                                          [(1, self.hplace[h][k]) for k in range(self.hrank(h, i)+1)],
                                          ">=", self.hosp_cap[h], "Hosp has space var has correct value"))
                var = self.pb_model.create_var("type1-{}-{}".format(i, j))
                self.bp_vars.append(var)
                self.pb_model.add_constr(Constraint([(1, var), (-1, hosp_has_space_var)] +
                                          [(1, v) for v in self.rplace[i][:j+1]], ">=", 0, "Type 1 stability"))

        # TYPE 2
        for i, j in self.couples:
            self.add_type2(i, j)
            self.add_type2(j, i)

        # TYPE 3
        for i, j in self.couples:
            self.add_type3a(i, j)
            self.add_type3bcd(i, j) 
            self.add_type3bcd(j, i)

        self.pb_model.add_sum_leq_constr(self.bp_vars, self.max_bp, "Max permitted number of blocking pairs")

    def add_type2(self, i, partner):
        for j, h in enumerate(self.rpref[i]):
            hrank_of_res = self.hrank(h, i)
            hrank_of_partner = self.hrank(h, partner) 
            hosp_has_space_var = self.pb_model.create_var("hosp_space-{}-{}".format(h, hrank_of_res))
            hplace_vars = [self.hplace[h][k] for k in range(hrank_of_res+1)]
            if hrank_of_partner is not None and hrank_of_partner > hrank_of_res:
                hplace_vars.append(self.hplace[h][hrank_of_partner])
            self.pb_model.add_constr(Constraint([(self.hosp_cap[h], hosp_has_space_var)] + 
                                      [(1, v) for v in hplace_vars],
                                      ">=", self.hosp_cap[h], "Hosp has space var has correct value"))
            var = self.pb_model.create_var("type2-{}-{}".format(i, j))
            self.bp_vars.append(var)
            self.pb_model.add_constr(Constraint([(1, var), (-1, hosp_has_space_var)] +
                              [(-1, self.rplace[i][idx]) for idx in range(j+1, len(self.rpref[i])) if self.rpref[partner][j]==self.rpref[partner][idx]],
                              ">=", -1, "Type 2 stability"))
        
    def add_type3a(self, i, partner):
        for j, (h, h2) in enumerate(zip(self.rpref[i], self.rpref[partner])):
            if h == h2: continue

            hrank_of_res = self.hrank(h, i)
            hrank_of_partner = self.hrank(h2, partner) 
            hosp1_has_space_var = self.pb_model.create_var("hosp1_space-{}-{}".format(h, hrank_of_res))
            hplace_vars = [self.hplace[h][k] for k in range(hrank_of_res+1)]
            self.pb_model.add_constr(Constraint([(self.hosp_cap[h], hosp1_has_space_var)] + 
                                      [(1, v) for v in hplace_vars],
                                      ">=", self.hosp_cap[h], "Hosp 1 has space var has correct value"))
            hosp2_has_space_var = self.pb_model.create_var("hosp2_space-{}-{}".format(h2, hrank_of_partner))
            hplace_vars = [self.hplace[h2][k] for k in range(hrank_of_partner+1)]
            self.pb_model.add_constr(Constraint([(self.hosp_cap[h2], hosp2_has_space_var)] + 
                                      [(1, v) for v in hplace_vars],
                                      ">=", self.hosp_cap[h2], "Hosp 2 has space var has correct value"))

            var = self.pb_model.create_var("type3a-{}-{}".format(i, j))
            self.bp_vars.append(var)
            self.pb_model.add_constr(Constraint([(1, var), (-1, hosp1_has_space_var), (-1, hosp2_has_space_var)] +
                          [(-1, self.rplace[i][idx]) for idx in range(j+1, len(self.rpref[i]))
                                if h != self.rpref[i][idx] and h2 != self.rpref[partner][idx]] +
                          [(-1, self.r_unassigned[i])],
                          ">=", -2, "Type 3a stability"))
        
    def add_type3bcd(self, i, partner):
        # The hospital should have no more than capacity-2 spaces used up to resident i,
        # and no more than capacity-1 spaces used up to resident partner
        for j, (h, h2) in enumerate(zip(self.rpref[i], self.rpref[partner])):
            if h != h2: continue

            hrank_of_res = self.hrank(h, i)
            hrank_of_partner = self.hrank(h, partner) 
            if hrank_of_res > hrank_of_partner: continue

            hosp_has_space_part_1_var = self.pb_model.create_var("hosp(3bcd)_space-{}-{}-part-1".format(h, hrank_of_res))
            hplace_vars = [self.hplace[h][k] for k in range(hrank_of_res+1)]
            self.pb_model.add_constr(Constraint([(self.hosp_cap[h]-1, hosp_has_space_part_1_var)] + 
                                      [(1, v) for v in hplace_vars],
                                      ">=", self.hosp_cap[h]-1, "Hosp (3bcd) has space var part 1 has correct value"))

            hosp_has_space_part_2_var = self.pb_model.create_var("hosp(3bcd)_space-{}-{}-part-2".format(h, hrank_of_partner))
            hplace_vars = [self.hplace[h][k] for k in range(hrank_of_partner+1)]
            self.pb_model.add_constr(Constraint([(self.hosp_cap[h], hosp_has_space_part_2_var)] + 
                                      [(1, v) for v in hplace_vars],
                                      ">=", self.hosp_cap[h], "Hosp (3bcd) has space var part 2 has correct value"))

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
        
    def is_single(r):
        "Is resident r single?"
        return r >= self.first_single

    def rrank(self, r, h):
        """What ranks does resident r give hospital h (as an array)?
           Note that a resident in a couple may rank a hospital more than once
        """
        return [i for i, hosp in enumerate(self.rpref[r]) if hosp == h]
            
    def hrank(self, h, r):
        "What rank does hospital h give resident r?"
        for i, res in enumerate(self.hpref[h]):
            if res == r:
                return i
        return None
            
def main(lines, max_bp, quiet):
    instance = Instance(lines, PBModel(), max_bp)
    instance.write(quiet)

if __name__=="__main__":
    max_bp = int(sys.argv[1])
    quiet = len(sys.argv) > 2 and 'q' in sys.argv[2]
    main([line.strip() for line in sys.stdin.readlines() if line.strip()], max_bp, quiet)
