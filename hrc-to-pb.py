import sys
from pprint import pprint as pp
from collections import Counter

class PBModel(object):
    def __init__(self):
        self.var_names = []
        self.constrs = []

    def create_var(self, name):
        self.var_names.append(name)
        return len(self.var_names) - 1

    def add_constr(self, constr):
        self.constrs.append(constr)

    def add_sum_leq_constr(self, terms, rhs):
        self.constrs.append(([(1, term) for term in terms], "<=", rhs))

    def add_sum_geq_constr(self, terms, rhs):
        self.constrs.append(([(1, term) for term in terms], ">=", rhs))

    def add_sum_eq_constr(self, terms, rhs):
        self.constrs.append(([(1, term) for term in terms], "=", rhs))

    def show_constrs(self):
        for c in self.constrs:
            print "c ", " ".join("{}*{}".format(i, self.var_names[j]) for i, j in c[0]), c[1], c[2]

    def add_objective(self, objective):
        self.objective = objective

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
            self.pb_model.add_sum_eq_constr(self.rplace[-1], 1)
            self.rplace.append(self.rplace[-1])  # Identical list for other member of couple
        for i in range(self.first_single, self.nres):
            self.rplace.append([self.pb_model.create_var("res{}-{}".format(i, j))
                                    for j in range(len(self.rpref[i]) + 1)])
            # This resident is assigned exactly one position on her pref list (or unassigned)
            self.pb_model.add_sum_eq_constr(self.rplace[-1], 1)

        # For each hosp, an array of vars such that
        # self.hplace[i][j]==1 <-> hospital i gets its j^th choice
        # (with the last position indicting not full)
        self.hplace = []
        for i in range(self.nhosp):
            self.hplace.append([self.pb_model.create_var("hosp{}-{}".format(i, j))
                                    for j in range(len(self.hpref[i]) + 1)])
            # Hospital capacity constraint
            self.pb_model.add_sum_leq_constr(self.hplace[-1], self.hosp_cap[i])
            # Force "unmatched" var if an insufficient number of places are matched
            self.pb_model.add_constr(([(1, j) for j in self.hplace[-1][:-1]] +
                                      [(self.hosp_cap[i], self.hplace[-1][-1])], ">=",
                                      self.hosp_cap[i]))

        # Chosen hosp prefs match chosen res prefs
        for i, prefs in enumerate(self.hpref):
            for pos, res in enumerate(prefs):
                self.pb_model.add_constr(([(1, self.hplace[i][pos])] +
                                 [(-1, self.rplace[res][k]) for k in self.rrank(res, i)], "=", 0))


        self.add_stability()
            
#        pp(self.rplace)
#        pp(self.hosp_cap)
#        pp(self.hplace)
        self.pb_model.show_constrs()

    def add_stability(self):
        self.bp_vars = []   # Blocking pair vars
        for i in range(self.first_single, self.nres):
            for j, h in enumerate(self.rpref[i]):
                hosp_has_space_var = self.pb_model.create_var("hosp_space-{}-{}".format(h, self.hrank(h, i)))
                self.pb_model.add_constr(([(1, hosp_has_space_var)] + 
                                          [(1, self.hplace[h][k]) for k in range(self.hrank(h, i))],
                                          "<=", self.hosp_cap[h]))
                var = self.pb_model.create_var("type1-{}-{}".format(i, j))
                self.pb_model.add_constr(([(1, var), (-1, hosp_has_space_var)] +
                                          [(1, v) for v in self.rplace[i][:j+1]], ">=", 0))

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
