"""
    Builds LAMMPS LJ Force Field files from a csv including individual atom parameters.
        This enforces mixing rules, & Drude Symmetrization Conditions

    Written By: Amro Dodin
"""
import pandas as pd
import mixing_rules as mr


default_db = ['./ff.csv']


class ForceField:
    def __init__(self, type_list, polarizable=False, symmetrized=True, db_fnames=default_db,
                 mixing_rules=(mr.arithmetic, mr.geometric, mr.geometric), drude_mass=0.8):
        # Collect Force Field Databases
        dbs = []
        for dbname in db_fnames:
            db = pd.read_csv(dbname, index_col=0)
            dbs.append(db)
        self._db = pd.concat(dbs)

        # Store Force Field Type Specification
        self._polarizable = polarizable
        self._mixing_rules = mixing_rules
        self._drude_mass = drude_mass
        self._symmetrized = symmetrized

        # Construct List of atom Types
        self._build_type_list(type_list)

        # Build Polarizable Force Field Details
        self._build_polarizable_types()
        self._num_types = len(self._type_list) + len(self._polarizable_types)

        # Build Force Field
        self._build_masses()
        self._build_coefficients()
        self._build_refs()

    def get_spec(self):
        if self._polarizable:
            return {'polarizable': self._polarizable,
                    'symmetrized': self._symmetrized,
                    'mixing_rules': self._mixing_rules,
                    'drude_mass': self._drude_mass}
        else:
            return {'polarizable': self._polarizable,
                    'mixing_rules': self._mixing_rules}

    def get_type_list(self):
        return self._type_list

    def __repr__(self):
        return 'Masses:\n' + str(self._masses) + '\nCoefficients:\n' + str(self._coefficients) + '\nReferences:\n' + str(self._refs)

    def _build_polarizable_types(self):
        self._polarizable_types = []
        self._pol_string = ''
        if self._polarizable:
            for t in self._type_list:
                if self._db.loc[t, 'alpha'] > 0:
                    self._polarizable_types.append(t)
                    self._pol_string += 'C '
                else:
                    self._pol_string += 'N '
            for t in self._polarizable_types:
                self._pol_string += 'D '

    def _mass_ratio(self, label):
        if label[-1] == 'S':
            return self._drude_mass / self._db.loc[label[:-1], 'mass']
        else:
            return 1 - self._drude_mass / self._db.loc[label[:-1], 'mass']

    def _build_type_list(self, at_mol_list):
        """ Builds a list of types from list of atom & molecule names

        :param at_mol_list: List of atom and Molecule names (e.g. ['SPC',  'Na', 'Cl' ])
        :return: list of atom types in order provided (e.g. ['OSPC', 'HSPC', 'Na', 'Cl'])
        """
        self._type_list = []
        for t in at_mol_list:
            if t in self._db.Molecule.unique():
                loc_list = self._db[self._db.Molecule == t].index.unique()
                for a in loc_list:
                    self._type_list.append(a)
            else:
                self._type_list.append(t)

    def _build_coefficients(self):
        tuples, labels = [], []
        sigmas, epsilons = [], []
        if self._polarizable_types:
            alphas = []
        for i in range(self._num_types):
            at_i, label_i = identify_pol_type(i, self._type_list, self._polarizable_types)
            if at_i in self._polarizable_types:
                mass_ratio_i = self._mass_ratio(label_i)
            else:
                mass_ratio_i = 1
            for j in range(i, self._num_types):
                at_j, label_j = identify_pol_type(j, self._type_list, self._polarizable_types)
                if at_j in self._polarizable_types:
                    mass_ratio_j = self._mass_ratio(label_j)
                else:
                    mass_ratio_j = 1
                tuples.append((i + 1, j + 1))
                labels.append("%s-%s" % (label_i, label_j))
                sigmas.append(self._mixing_rules[0](self._db.loc[at_i, 'sigma'], self._db.loc[at_j, 'sigma']))
                if self._symmetrized:
                    epsilons.append(
                        mr.drude_symmetrized(self._db.loc[at_i, 'epsilon'], self._db.loc[at_j, 'epsilon'], mass_ratio_i, mass_ratio_j,
                                             self._mixing_rules[1])
                    )
                else:
                    epsilons.append(
                        mr.drude_asymmetric(self._db.loc[at_i, 'epsilon'], self._db.loc[at_j, 'epsilon'], mass_ratio_i,
                                             mass_ratio_j,
                                             self._mixing_rules[1])
                    )
                if self._polarizable_types:
                    alphas.append(self._mixing_rules[2](self._db.loc[at_i, 'alpha'], self._db.loc[at_j, 'alpha']))
        if self._polarizable_types:
            data = {'sigma': sigmas, 'epsilon': epsilons, 'alpha': alphas, 'label': labels}
        else:
            data = {'sigma': sigmas, 'epsilon': epsilons, 'label': labels}
        index = pd.MultiIndex.from_tuples(tuples)
        self._coefficients = pd.DataFrame(data, index=index)

    def _build_masses(self):
        labels = []
        masses = []
        for t in self._type_list:
            if t not in self._polarizable_types:
                labels.append(t)
                masses.append(self._db.loc[t, 'mass'])
            else:
                labels.append(t + 'C')
                masses.append(self._db.loc[t, 'mass'] - self._drude_mass)
        for t in self._polarizable_types:
            labels.append(t + 'D')
            masses.append(self._drude_mass)
        self._masses = pd.DataFrame({'mass': masses, 'label': labels})
        self._masses.index += 1

    def _build_refs(self):
        self._refs = self._db.loc[self._type_list, 'Reference']

    def print(self, fname=None):
        if fname:
            f = open(fname, 'w')
            printer = lambda x: f.write(x + '\n')
        else:
            printer = print

        # Print Header
        printer('# LAMMPS Force Field File')
        printer('# Generated by ForceField.py')
        printer('')

        # Print Reference List
        printer('# FF References:')
        for i, r in self._refs.iteritems():
            printer('#      %s: %s' % (i, r))
        printer('')

        # Print The Drude Fix line if polarizable
        if self._pol_string:
            printer('fix drude all drude ' + self._pol_string)
            printer('pair_style hybrid/overlay lj/cut/coul/long 12.0 thole 2.6 12.0')
            pair_style = 'lj/cut/coul/long'
        else:
            printer('pair_style lj/cut/coul/long 12.0')
            pair_style = ''
        printer('kspace_style pppm 1.0e-3')
        printer('')

        # Print Masses
        printer('# Masses:')
        for i, r in self._masses.iterrows():
            printer('mass %2d %5.3f # %s' % (i, r.mass, r.label))
        printer('')

        # Print FF Coefficients
        printer('# Pair Potential Coefficients')
        for i, r in self._coefficients.iterrows():
            printer('pair_coeff %2d %2d %s %8.5f %8.5f # %s' % (i[0], i[1], pair_style, r.epsilon, r.sigma, r.label))
        printer('')

        # Print Thole Screening
        if self._pol_string:
            printer('# Thole Screening')
            for i, r in self._coefficients.iterrows():
                if r.alpha > 0:
                    printer('pair_coeff %2d %2d thole %8.5f # %s' % (i[0], i[1], r.alpha, r.label))

        if fname:
            f.close()


# Utility Functions
def identify_pol_type(i, type_list, polarizable_types):
    """ Identifies which atom types are drude Cores vs Shells & generates labels

    :param i: index of type to check
    :param type_list: list of physical atom types
    :param polarizable_types: list of polarizable atom types
    :return: atom name and label (e.g. O -> OC or OS)
    """
    num_elements = len(type_list)
    if i < num_elements:
        ati = type_list[i]
        if ati in polarizable_types:
            labeli = ati + 'C'
        else:
            labeli = ati
    else:
        ati = polarizable_types[num_elements - i]
        labeli = ati + 'S'
    return ati, labeli

