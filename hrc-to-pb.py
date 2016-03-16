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

    def show_constrs(self):
        for c in self.constrs:
            print "* ", c.name+":"
            print "*     ", " ".join("{}*{}".format(i, self.var_names[j]) for i, j in c[0]), c[1], c[2]

    def show_objective(self):
        print "* Objective: max:"
        print "*     ", " ".join("{}*{}".format(i, self.var_names[j]) for i, j in self.objective)

    def show_var_names(self):
        for i, name in enumerate(self.var_names):
            print "*", i+1, name

    def add_objective(self, objective):
        self.objective = objective

    def write_model_size_comment(self):
        print "* #variable= {} #constraint= {}".format(len(self.var_names), len(self.constrs))

    def write_model(self):
        print "min:", " ".join("-" + str(t[0]) + " " + "x{}".format(t[1]+1) for t in self.objective) + ";"
        for c in self.constrs:
            print "* ", c.name+":"
            print "*     ", " ".join("{}*{}".format(i, self.var_names[j]) for i, j in c[0]), c[1], c[2]
            if c.comp == "<=":
                print " ".join("{:+d} x{}".format(-i, j+1) for i, j in c[0]), ">=", str(c[2]) + ";"
            else:
                print " ".join("{:+d} x{}".format(i, j+1) for i, j in c[0]), c[1], str(c[2]) + ";"
            print "*"
            

class Instance(object):
    def __init__(self, lines, pb_model):
        self.pb_model = pb_model

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
#            # Hospital capacity constraint
#            self.pb_model.add_sum_leq_constr(self.hplace[-1], self.hosp_cap[i], "Hospital capacity")
#            # Force "unmatched" var if an insufficient number of places are matched
#            self.pb_model.add_constr(Constraint([(1, j) for j in self.hplace[-1][:-1]] +
#                                                [(self.hosp_cap[i], self.hplace[-1][-1])], ">=",
#                                                self.hosp_cap[i], "Hosp unmatched"))

#        # Chosen hosp prefs match chosen res prefs
#        for i, prefs in enumerate(self.hpref):
#            for pos, res in enumerate(prefs):
#                self.pb_model.add_constr(
#                        Constraint([(1, self.hplace[i][pos])] +
#                                   [(-1, self.rplace[res][k]) for k in self.rrank(res, i)],
#                                   "=", 0, "Hosp pref matches res prefs"))


#        self.add_stability()
            
        self.add_objective()
#        pp(self.rplace)
#        pp(self.hosp_cap)
#        pp(self.hplace)
        self.pb_model.write_model_size_comment()
#        self.pb_model.show_constrs()
        self.pb_model.show_var_names()
        self.pb_model.show_objective()
        self.pb_model.write_model()

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
        for i in range(self.first_single, self.nres):
            for j, h in enumerate(self.rpref[i]):
                hosp_has_space_var = self.pb_model.create_var("hosp_space-{}-{}".format(h, self.hrank(h, i)))
                self.pb_model.add_constr(Constraint([(1, hosp_has_space_var)] + 
                                          [(1, self.hplace[h][k]) for k in range(self.hrank(h, i))],
                                          "<=", self.hosp_cap[h], "Hosp has space var has correct value"))
                var = self.pb_model.create_var("type1-{}-{}".format(i, j))
                self.bp_vars.append(var)
                self.pb_model.add_constr(Constraint([(1, var), (-1, hosp_has_space_var)] +
                                          [(1, v) for v in self.rplace[i][:j+1]], ">=", 0, "Type 1 stability"))

    def rrank(self, r, h):
        """What ranks does resident r give hospital h (as an array)?
           Note that a resident in a couple may rank a hospital more than once
        """
        ranks = []
        for i, hosp in enumerate(self.rpref[r]):
            if hosp == h:
                ranks.append(i)
        if ranks:
            return ranks
        else:
            return [len(rpref[r])]
            
    def hrank(self, h, r):
        "What rank does hospital h give resident r?"
        for i, res in enumerate(self.hpref[h]):
            if res == r:
                return i
        return len(hpref[h])
            


def main(lines):
    instance = Instance(lines, PBModel())


if __name__=="__main__":
    main([line.strip() for line in sys.stdin.readlines() if line.strip()])
