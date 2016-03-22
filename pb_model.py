from collections import namedtuple

Constraint = namedtuple('Constraint', ['terms', 'comp', 'rhs', 'name'])

class PBModel(object):
    def __init__(self, flatzinc):
        self.var_names = []
        self.constrs = []
        self.flatzinc = flatzinc

    def create_var(self, name):
        self.var_names.append(name)
        return len(self.var_names) - 1

    def add_constr(self, constr):
        self.constrs.append(constr)

    def add_sum_leq_constr(self, variables, rhs, name="UNNAMED"):
        self.constrs.append(Constraint([(1, variable) for variable in variables], "<=", rhs, name))

    def add_sum_eq_constr(self, variables, rhs, name="UNNAMED"):
        self.constrs.append(Constraint([(1, variable) for variable in variables], "=", rhs, name))

    def add_exactly_one_constr(self, variables, name="UNNAMED"):
        self.constrs.append(Constraint([(1, variable) for variable in variables], "=", 1, name))

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
            
    def write_flatzinc(self):
        print "array[1..{}] of var 0..1: x;".format(len(self.var_names))
        print "var int: obj :: output_var;"
        print "constraint int_lin_eq([{},-1],[{},obj],0);".format(
                ",".join(str(t[0]) for t in self.objective),
                ",".join("x[{}]".format(t[1]+1) for t in self.objective))
        for c in self.constrs:
            if c.comp == "<=":
                print "constraint int_lin_le([{}],[{}],{});".format(
                        ",".join(str(t[0]) for t in c.terms),
                        ",".join("x[{}]".format(t[1]+1) for t in c.terms),
                        c.rhs)
            elif c.comp == ">=":
                print "constraint int_lin_le([{}],[{}],{});".format(
                        ",".join(str(-t[0]) for t in c.terms),
                        ",".join("x[{}]".format(t[1]+1) for t in c.terms),
                        -c.rhs)
            elif c.comp == "=":
                print "constraint int_lin_eq([{}],[{}],{});".format(
                        ",".join(str(t[0]) for t in c.terms),
                        ",".join("x[{}]".format(t[1]+1) for t in c.terms),
                        c.rhs)
        print "solve maximize obj;"
            
    def write(self, quiet):
        if self.flatzinc:
            self.write_flatzinc()
        else:  # .opb
            self.write_model_size_comment()
            if not quiet:
                self.show_var_names()
                self.show_objective()
            self.write_model(quiet)

