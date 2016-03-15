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
            self.rplace.append([self.pb_model.create_var("res{}-{}".format(i, j))
                                    for j in range(len(self.rpref[i]) + 1)])
            self.rplace.append(self.rplace[-1])  # Identical list for other member of couple
        for i in range(self.first_single, self.nres):
            self.rplace.append([self.pb_model.create_var("res{}-{}".format(i, j))
                                    for j in range(len(self.rpref[i]) + 1)])

        # For each hosp, an array of vars such that
        # self.hplace[i][j]==1 <-> hospital i gets its j^th choice
        # (with the last position indicting not full)
        self.hplace = []
        for i in range(self.nhosp):
            self.hplace.append([self.pb_model.create_var("hosp{}-{}".format(i, j))
                                    for j in range(len(self.hpref[i]) + 1)])
            
#        pp(self.rplace)
#        pp(self.hosp_cap)
#        pp(self.hplace)

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
