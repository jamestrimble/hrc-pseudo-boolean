import sys
from pprint import pprint as pp
from collections import Counter, namedtuple

Constraint = namedtuple('Constraint', ['terms', 'comp', 'rhs', 'name'])

class PBModel(object):
    def __init__(self):
        self.var_names = []
        self.constrs = []

    def create_var(self, name):
        self.var_names.append(name)
        return len(self.var_names) - 1

    def add_constr(self, constr):
        self.constrs.append(constr)

    def add_sum_leq_constr(self, terms, rhs, name="UNNAMED"):
        self.constrs.append(Constraint([(1, term) for term in terms], "<=", rhs, name))

    def add_sum_geq_constr(self, terms, rhs, name="UNNAMED"):
        self.constrs.append(Constraint([(1, term) for term in terms], ">=", rhs, name))

    def add_sum_eq_constr(self, terms, rhs, name="UNNAMED"):
        self.constrs.append(Constraint([(1, term) for term in terms], "=", rhs, name))

    def show_objective(self):
        print "* Objective: max:"
        print "*     ", " ".join("{}*{}".format(i, self.var_names[j]) for i, j in self.objective)

    def show_var_names(self):
        for i, name in enumerate(self.var_names):
            print "*", "NAME", i+1, name

    def add_objective(self, objective):
        self.objective = objective

    def write_model_size_comment(self):
        print "* #variable= {} #constraint= {}".format(len(self.var_names), len(self.constrs))

    def write_model(self, quiet):
        print "min:", " ".join("-" + str(t[0]) + " " + "x{}".format(t[1]+1) for t in self.objective) + ";"
        for c in self.constrs:
            if not quiet:
                print "* ", c.name+":"
                print "*     ", " ".join("{}*{}".format(i, self.var_names[j]) for i, j in c[0]), c[1], c[2]
            if c.comp == "<=":
                print " ".join("{:+d} x{}".format(-i, j+1) for i, j in c[0]), ">=", str(-c[2]) + ";"
            else:
                print " ".join("{:+d} x{}".format(i, j+1) for i, j in c[0]), c[1], str(c[2]) + ";"
            if not quiet: print "*"
            

class Instance(object):
    def __init__(self, lines, pb_model, max_bp):
        self.pb_model = pb_model
        self.max_bp = max_bp   # Maximum permitted number of blocking pairs

        self.nres = int(lines[0])
        self.nhosp = int(lines[1])
        self.ncoup = int(lines[2])
        self.nsingle = self.nres - 2 * self.ncoup
        self.first_single = 2 * self.ncoup
        self.npost = int(lines[3])
        lines = lines[9:]

        # Pref list for each resident
        self.rpref = [[int(x) for x in line.split()[1:]] for line in lines[:self.nres]]
        
        # Pref list for each hospital
        self.hpref = [[int(x) for x in line.split()[2:]] for line in lines[self.nres:self.nres+self.nhosp]]

        # Hospital capacities
        self.hosp_cap = [int(line.split()[1]) for line in lines[self.nres:self.nres+self.nhosp]]

        # For each resident, an array of vars such that
        # self.rplace[i][j]==1 <-> resident i gets his j^th choice
        # (with the last position indicting unplaced)
        self.rplace = []
        for i in range(self.ncoup):
            self.rplace.append([self.pb_model.create_var("res{}-{}".format(i*2, j))
                                    for j in range(len(self.rpref[i*2]) + 1)])
            # This resident is assigned exactly one position on her pref list (or unassigned)
            self.pb_model.add_sum_eq_constr(self.rplace[-1], 1, "Coupled resident assigned only once")
            self.rplace.append(self.rplace[-1])  # Identical list for other member of couple
        for i in range(self.first_single, self.nres):
            self.rplace.append([self.pb_model.create_var("res{}-{}".format(i, j))
                                    for j in range(len(self.rpref[i]) + 1)])
            # This resident is assigned exactly one position on her pref list (or unassigned)
            self.pb_model.add_sum_eq_constr(self.rplace[-1], 1, "Single resident assigned only once")

        # For each hosp, an array of vars such that
        # self.hplace[i][j]==1 <-> hospital i gets its j^th choice
        # (with the last position indicting not full)
        self.hplace = []
        for i in range(self.nhosp):
            self.hplace.append([self.pb_model.create_var("hosp{}-{}".format(i, j))
                                    for j in range(len(self.hpref[i]) + 1)])
            # Hospital capacity constraint
            self.pb_model.add_sum_leq_constr(self.hplace[-1], self.hosp_cap[i], "Hospital capacity")
            # Force "unmatched" var if an insufficient number of places are matched
            self.pb_model.add_constr(Constraint([(1, j) for j in self.hplace[-1][:-1]] +
                                                [(self.hosp_cap[i], self.hplace[-1][-1])], ">=",
                                                self.hosp_cap[i], "Hosp unmatched"))

        # Chosen hosp prefs match chosen res prefs
        for i, prefs in enumerate(self.hpref):
            for pos, res in enumerate(prefs):
                self.pb_model.add_constr(
                        Constraint([(1, self.hplace[i][pos])] +
                                   [(-1, self.rplace[res][k]) for k in self.rrank(res, i)],
                                   "=", 0, "Hosp pref matches res prefs"))


        self.add_stability()
            
        self.add_objective()

    def write(self, quiet):
        self.pb_model.write_model_size_comment()
        if not quiet:
            self.pb_model.show_var_names()
        self.pb_model.show_objective()
        self.pb_model.write_model(quiet)

    def add_objective(self):
        obj_terms = []
        for i in range(self.ncoup):
            for v in self.rplace[i*2][:-1]:
                obj_terms.append((2, v))
        for i in range(self.nsingle):
            for v in self.rplace[self.first_single + i][:-1]:
                obj_terms.append((1, v))
        self.pb_model.add_objective(obj_terms)
        
    def add_stability(self):
        self.bp_vars = []   # Blocking pair vars

        # TYPE 1
        for i in range(self.first_single, self.nres):
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
        for couple_num in range(self.ncoup):
            i = couple_num * 2
            self.add_type2(i, i+1)
            self.add_type2(i+1, i)

        # TYPE 3
        for couple_num in range(self.ncoup):
            i = couple_num * 2
            self.add_type3a(i, i+1)
            self.add_type3bcd(i, i+1)
            self.add_type3bcd(i+1, i)

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
                              [(-1, self.rplace[i][idx]) for idx in range(j+1, len(self.rpref[i])+1)
                                    if idx == len(self.rpref[i]) or (h != self.rpref[i][idx] and h2 != self.rpref[partner][idx])],
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
            self.pb_model.add_constr(Constraint([(1, var), (-1, hosp_has_space_var)] +
                              [(-1, self.rplace[i][idx]) for idx in range(j+1, len(self.rpref[i])+1)
                                    if idx == len(self.rpref[i]) or (h != self.rpref[i][idx] and h != self.rpref[partner][idx])],
                              ">=", -1, "Type 3bcd stability"))
        
    def rrank(self, r, h):
        """What ranks does resident r give hospital h (as an array)?
           Note that a resident in a couple may rank a hospital more than once
        """
        ranks = []
        for i, hosp in enumerate(self.rpref[r]):
            if hosp == h:
                ranks.append(i)
        return ranks
            
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
