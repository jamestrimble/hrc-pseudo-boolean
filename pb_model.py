from collections import namedtuple

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
            

