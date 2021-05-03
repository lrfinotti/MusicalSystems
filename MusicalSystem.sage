'''
Provides functions.  FIXME!
'''

from itertools import combinations


def vec_to_string(v, delim='[]'):
    '''
    Converts a given vector to a string suitable for printing.
    Assumes entries are up to two digits for correct spacing.

    The optional argument 'delim' gives option of delimiters.

    This is used to print tables.

    EXAMPLE:
    sage: vec_to_string([1,2,3,4])
    '[ 1,  2,  3,  4]'

    sage: vec_to_string([1,2,3,4], delim='()')
    '( 1,  2,  3,  4)'
    '''
    vx = [f'{x:>2}' for x in v]

    return delim[0] + ' ,'.join(vx) + delim[-1]


def print_prime_form(v, letters='te'):
    '''
    Converts a given vector to a string in prime form format.

    Giving letters as 'TE' makes it use capital letters.


    EXAMPLE:
    sage: vec_to_string([1,2,10,11])
    '(12te)'

    sage: vec_to_string([1,2,10,11], letters='TE')
    '(12te)'
    '''
    res = '('

    for x in v:
        if 0 <= x <= 9:
            res += str(x)
        elif x == 10:
            res += letters[0]
        elif x == 11:
            res += letters[1]

    res += ')'

    return res


class MusicalSystem:
    '''
    Given two permutations, rho of order 12, and phi of order 2,  generating
    a Dihedral group, creates the corresponding "musical system".  By default,
       rho = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)   and
       phi = (1, 11), (2, 10), (3, 9), (4, 8), (5, 7),
    which gives the standard system.
    '''

    def __init__(self, rho=None, phi=None):
        S12 = SymmetricGroup(range(12))

        # define rho and phi
        if rho:
            self.rho = rho
        else:
            self.rho = S12([(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)])

        if phi:
            # try to find self.phi fixing 0 (always works if transitive!)
            tphi = phi
            for i in range(12):
                if tphi(0) == 0:
                    self.phi = tphi
                    break
                tphi *= self.rho
            else: # not found!
                self.phi = phi
        else:
            self.phi = S12([(1, 11), (2, 10), (3, 9), (4, 8), (5, 7)])

        self.group = S12.subgroup([self.rho, self.phi])

        # check if Dihedral
        if not (self.group).is_isomorphic(DihedralGroup(12)):
            raise ValueError(
                'Given permuations do not give Dihedral Group of order 24.')

        # is rho a 12-cycle?
        self.is_transitive = (len(self.rho.orbit(0)) == 12)

        # find sort/interval vector
        tmpv = list(range(12))
        svec = []
        while tmpv:
            x = tmpv[0]
            while x not in svec:
                svec.append(x)
                tmpv.remove(x)
                x = self.rho(x)

        self.intvec = svec

    @classmethod
    def from_conj(cls, sigma):
        '''
        Creates a MusicalSystem from conjugation.  Given a permutaiton
        sigma is S12, creates the system with group sigma*D12*sigma^(-1).
        '''
        S12 = SymmetricGroup(range(12))
        rho = S12([(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)])
        phi = S12([(1, 11), (2, 10), (3, 9), (4, 8), (5, 7)])
        rhoh = sigma ^ (-1)*rho*sigma
        phih = sigma ^ (-1)*phi*sigma
        return cls(rhoh, phih)

    def __repr__(self):
        return f'MusicalSystem(rho={self.rho}, phi={self.phi}): Musical System given by rotation {self.rho} and reflection {self.phi}'

    def __str__(self):
        return f'MusicalSystem(rho={self.rho}, phi={self.phi})'

    def __eq__(self, other):
        if not isinstance(other, MusicalSystem):
            return False
        # rho is important to measure intervals!
        if self.rho != other.rho:
            return False
        # does phi matter, or just the group?
        # if self.phi != other.phi:
        #     return False
        return self.group == other.group

    def absolute_interval(self, a, b):
        '''
        Gives the ABSOLUTE interval, measure by rho, using the
        self.intvec
        '''
        svec = self.intvec
        return (svec.index(b) - svec.index(a)) % 12

    def interval(self, a, b):
        '''
        Gives the interval invariant by rho!  (Important when rho is not
        12 cycle!)
        '''
        if self.is_transitive:
            return self.absolute_interval(a, b)
        else:
            rho = self.rho
            return min((self.absolute_interval((rho ^ i)(a), (rho ^ i)(b)) for i in range(12)))

    def ints(self, v):
        '''
        Given a vector v, gives a a list of intervals v[0] to v[-1], v[0] to v[-2],
        etc.  Useful for prime and normal forms!
        '''
        return [self.interval(v[0], v[-i]) for i in range(1, len(v))]

    def interval_class(self, a, b):
        '''
        Gives the interval invariant under the action of the group.
        '''
        phi = self.phi
        return min(self.interval(a, b), self.interval(phi(a), phi(b)))

    def sveckey(self, x):
        '''
        Gives the key function to sort according to self.intvec
        '''
        return (self.intvec).index(x)

    def first(self, v, w):
        '''
        Given two vectors, see which is first with the ordering given
        by self.intvec.
        '''
        if v == w:
            return v

        for x, y in zip(v, w):
            if self.sveckey(x) < self.sveckey(y):
                return v
            if self.sveckey(x) > self.sveckey(y):
                return w
            return None  # something went worng if here

    def sort_vec(self, v):
        '''
        Sort a vector according to self.intvec.

        First, it start with 0.  For the next element, we use the
        one which needs the least applitcations of rho to get to the
        next.

        If rhoh = rho = (0,1,2...,11), then this is the regular sort.
        '''
        svec = self.intvec
        res = []

        for x in svec:
            if x in v:
                res.append(x)
        return res

    def normal_form(self, v):
        '''
        Orders the entries of v so that they are most compressed
        to the left.  (Not necessarily starting at the first)
        '''
        w = self.sort_vec(v)
        res = w
        l = len(v)

        for _ in range(1, l):
            w = w[1:] + [w[0]]
            lres = self.ints(res)
            lw = self.ints(w)
            if lw < lres:
                res = w

        return res

    # SEEMS TO NOT BE USED!
    def min_vec(self, v, w):
        '''
        Given two vectors, finds which one comes first, i.e., it
        is more "concentrated to the left".  This means that the
        last element of the one that comes first is obtained from
        0 with less applications of rhoh than the one that comes after.
        (If rhoh = rho, it means that the last element of the vector
        that comes first is less than the last element of the vector
        that comes second.)

        If the last matches, proceeds to one before, and so on.

        ASSUMES: v and w are sorted, v[0] = w[0], and same length.

        '''
        if v == w:
            return v
        vv = list(v)  # copy of v
        ww = list(w)  # copy of w

        if vv[0] != ww[0]:
            raise ValueError('The first coordinates must be equal.')

        if self.ints(v) < self.ints(w):
            return v
        else:
            return w

    def symmetries(self, v, w=None):
        '''
        Gives the powers of rho and phi that map v to w as sets.
        If only v is given, check inner symmetries.  (Always has
        identity as one.)
        '''
        if w is None:
            w = list(v)
        else:
            w = list(w)  # do not change original!

        u = list(v)  # changed v (copy)

        if len(v) != len(w):
            return []

        res = []
        sw = set(w)

        for i in range(12):
            if set(u) == sw:
                res.append((i, 0))
            u = self.vrho(u)

        u = self.vphi(v)
        for i in range(12):
            if set(u) == sw:
                res.append((i, 1))
            u = self.vrho(u)

        return res

    def symmetry_maps(self, v, w):
        '''
        Lists all symmetries that map v to w.  The result is given in pairs
        with the first element the power of rho, and the second the power of
        phi.
        '''
        if len(v) != len(w):
            return []

        res = set()
        for r in range(12):
            for s in range(2):
                if self.vmap(r, s, v, sort='sort') == sorted(w):
                    res.add((r,s))
        return list(res)


    def class_sum(self, v):
        '''
        Returns the class sum of the pitch class set v.
        '''
        return sum(self.interval(0, x) for x in v) % 12


    def find_first_webern_rows(self, v, inner_sym=True):
        '''
        Given a set class v of 3 elements, finds other 3 pitch classes in the
        same set class that covers all 12 pitch classes.

        It can be used as first row of the series, to produce a series that
        preserve the set classes in sets of 3, like Webern.

        This can still be tweaked by changing order in the sets of 3, the order
        of the sets of 3, applying rhoh^i to all, etc.
        '''
        if not self.is_transitive:  # not a 12-cycle!
            raise TypeError("rho is not a 12-cycle.")

        rho = self.rho
        phi = self.phi

        vv = set()  # all elements in the set class
        w = tuple(sorted(v))
        for _ in range(0, 12):
            vv.add(w)
            # w = tuple(sorted([rho(x) for x in w]))
            w = tuple(self.vrho(w, sort='sort'))

        # w = tuple(sorted([phi(x) for x in v]))
        w = tuple(self.vphi(w, sort='sort'))
        for _ in range(0, 12):
            vv.add(w)
            # w = tuple(sorted([rho(x) for x in w]))
            w = tuple(self.vrho(w, sort='sort'))

        res = []
        vv = list(vv)

        # brute force to find all that work
        for X, Y, Z, T in combinations(vv, 4):
            if len(set(X).union(set(Y), set(Z), set(T))) == 12:
                if not inner_sym:
                    res.append(sorted([X, Y, Z, T]))
                    res.sort()
                else:
                    for W in [Y, Z, T]:
                        W1, W2 = [WW for WW in [Y, Z, T] if WW != W]
                        sc = SetClass(X + W, self)
                        nsym = sc.nsym()
                        if nsym[0] > 0 and nsym[1] > 0 and SetClass(W1 + W2, self) == sc:
                            res.append([X, W, W1, W2, nsym])
                    res.sort()
        return res

    def find_all_first_webern_rows(self, inner_sym=True):
        '''
        Find all possible first rows for a Webern matrix.
        '''
        if not self.is_transitive:  # not a 12-cycle!
            raise TypeError("rho is not a 12-cycle.")

        class_elements = ClassElements(3, self).matrix
        set_classes = [x[0] for x in class_elements]

        res = []
        for v in set_classes:
            res += self.find_first_webern_rows(v, inner_sym=inner_sym)
        return res

    def vrho(self, v, power=1, sort='normal'):
        '''
        Apply rho to the entries of a vector then take normal form.
        '''
        if sort == 'sort':
            return sorted(list(map(self.rho^power, v)))
        if sort == 'rho':
            return self.sort_vec(list(map(self.rho^power, v)))
        if sort == 'none':
            return list(map(self.rho^power, v))
        # if sort == 'normal' or anything else
        return self.normal_form(list(map(self.rho^power, v)))

    def vphi(self, v, power=1, sort='normal'):
        '''
        Apply phi to the entries of a vector then take normal form.
        '''
        if power % 2 == 0:
            if sort == 'sort':
                return sorted(v)
            if sort == 'rho':
                return self.sort_vec(v)
            if sort == 'none':
                return v
            # sort == 'normal' or anything else
            return self.normal_form(v)
        else:
            if sort == 'sort':
                return sorted(list(map(self.phi, v)))
            if sort == 'rho':
                return self.sort_vec(list(map(self.phi, v)))
            if sort == 'none':
                return list(map(self.phi, v))
            # if sort == 'normal' or anything else
            return self.normal_form(list(map(self.phi, v)))


    def vmap(self, r, s, v, sort='normal'):
        '''
        Apply rho^r*phi^s to entries of v.
        '''
        w = self.vphi(v, power=s, sort=sort)

        return self.vrho(w, power=r, sort=sort)


    def apply_map(self, map, v):
        '''
        Given an element of S_12, apply to all entries of v.
        '''
        if map not in S12:
            raise ValueError('Map must be in S_12')

        return [ map(x) for x in v ]

    def set_class(self, v):
        '''
        Creates the set class given by the pitch class set v for
        the system.
        '''
        return SetClass(v, MS=self)

    def forte_table(self, n):
        '''
        Creates the Forte Table for n-chords in the system.
        '''
        return ForteTable(n, MS=self)

    def webern_matrix(self, v0):
        '''
        Crteates the WebernMatrix with first row v0 in the system.
        '''
        return WebernMatrix(v0, MS=self)



class WebernMatrix:
    '''
    Given rhoh, phih and v0, gives the matrix, analog to the one from Webern,
    with labels on left, top, right and bottom.

    The default v0 is from Webern, so for arbitrary rhoh and phih, it will
    not preserve set classes.
    '''

    def __init__(self, v0=(0, 11, 3, 4, 8, 7, 9, 5, 6, 1, 2, 10), MS=MusicalSystem()):

        rho = MS.rho
        phi = MS.phi

        if not MS.is_transitive:  # not a 12-cycle!
            raise TypeError("rho is not a 12-cycle.")

        self.v0 = list(v0)
        self.rho = rho
        self.phi = phi
        self.MS = MS

        self.matrix = [v0] + [None for i in range(11)]
        self.left = [0] + [None for i in range(11)]  # left labels
        self.top = [None for i in range(12)]  # top labels

        w0 = [phi(x) for x in v0]

        # find power of rho^i*phi(v0) that starts with v0[0]
        # this will be first column
        # (gets top)
        for i in range(12):
            if w0[0] == v0[0]:
                self.top[0] = i
                break
            w0 = [rho(x) for x in w0]
        # else:
        #     raise NameError('Could not find w0 starting as v0.  Is rho a 12-cycle?')

        # find matrix with first entries matching the column w0 above
        # (gets self.left and matrix)
        for i in range(1, 12):
            v = [(rho ^ i)(x) for x in v0]
            j = w0.index(v[0])
            self.left[j] = i
            self.matrix[j] = v

        # find the powers of rho^j*phi that match the columns
        # (gets top)
        for i in (x for x in range(12) if x != self.top[0]):
            # remember Sage inverts operations
            w = [(phi*rho ^ i)(x) for x in v0]
            j = v0.index(w[0])
            self.top[j] = i

        # now we know that the last row ends with v0[0] and we can get the right
        self.right = [(self.left[i] - self.left[11]) % 12 for i in range(12)]

        # similarly we can get the bottom too
        self.bottom = [(self.top[i] + self.left[11]) % 12 for i in range(12)]

        vr = self.matrix[11]

        # ################## CHECK! ################## #
        # Verifies if everyting is working.
        # vr = self.matrix[11]
        # for i, j in enumerate(self.left):
        #     if self.matrix[i] != [(rho ^ j)(x) for x in v0]:
        #         raise NameError("Fail left!")

        # for i, j in enumerate(self.top):
        #     if [x[i] for x in self.matrix] != [(phi*rho ^ j)(x) for x in v0]:
        #         raise NameError("Fail top!")

        # for i, j in enumerate(self.right):
        #     if self.matrix[i] != [(rho ^ j)(x) for x in vr]:
        #         raise NameError("Fail right!")

        # for i, j in enumerate(self.bottom):
        #     if [x[i] for x in self.matrix] != [(phi*rho ^ j)(x) for x in vr]:
        #         raise NameError("Fail bottom")


    def __repr__(self):
        return 'Webern Table with first row {} and rho = {}, phi = {}'.format(self.v0, self.rho, self.phi)

    def __str__(self):
        S, left, top, right, bottom = self.matrix, self.left, self.top, self.right, self.bottom

        l = len(S)

        # top row
        vstr = ['  ', '|'] + ['{:>2}'.format(x) for x in top] + ['|']
        print('  '.join(vstr))
        mstr = '----|' + '--'*2*l + '--|----'
        # print(mstr)
        res = mstr+'\n'

        # middle rows
        for i, v in enumerate(S):
            vstr = ['{:>2}'.format(left[i]), '|'] + ['{:>2}'.format(x)
                                                     for x in v] + ['|', '{:>2}'.format(right[i])]
            # print('  '.join(vstr))
            res += '  '.join(vstr) + '\n'

        # bottom row
        # print(mstr)
        res += mstr + '\n'
        vstr = ['  ', '|'] + ['{:>2}'.format(x) for x in bottom] + ['|']
        # print('  '.join(vstr))
        res += '  '.join(vstr)

        return res

    def latex(self, srho='\\rho', sphi='\\phi', sid='1'):
        '''
        Print a LaTeX formatted table.   Use srho and sphi to name the functions, and
        sid for the identity/1.
        '''

        def fix_power(i, bphi=False):
            r'''
            Fixes the poowers when printing.  sid is the string for the identity.
            bphi tells if it is \rho^i \phi or just \rho^i.
            '''
            if i == 0:
                if bphi:
                    return f'${sphi}$'
                else:
                    return f'${sid}$'

            if i == 1:
                if bphi:
                    return f'${srho} {sphi}$'
                else:
                    return f'${srho}$'

            # i >= 2
            if bphi:
                return f'${srho}^{{{i}}} {sphi}$'
            else:
                return f'${srho}^{{{i}}}$'

        S, left, top, right, bottom = self.matrix, self.left, self.top, self.right, self.bottom

        res = '\\begin{tabular}{c|cccccccccccc|c}\n'
        # top row
        #vstr = ['  '] + [ '$\\hat{{\\rho}}^{{{}}}\\phi$'.format(x) for x in top  ] + [ '  ' ]
        vstr = ['  '] + [fix_power(x, bphi=True) for x in top] + ['  ']
        res += '  ' + ' &  '.join(vstr)
        res += '\\\\\n'
        res += '  \\hline\n'

        # middle rows
        for i, v in enumerate(S):
            vstr = [fix_power(left[i])] + \
                [f'${x}$' for x in v] + [fix_power(right[i])]
            res += '  ' + ' & '.join(vstr) + ' \\\\\n'

        res += '  \\hline\n'

        # bottom row
        vstr = ['  '] + [fix_power(x, bphi=True) for x in bottom] + ['  ']
        res += '  ' + ' &  '.join(vstr)
        res += '\n'
        res += '\\end{tabular}'

        return res

    def row(self, i):
        '''
        Gives the i-th row of the Webern matrix.
        '''
        return self.matrix[i]

    def column(self, j):
        '''
        Gives the j-th column of the Webern matrix:
        '''
        return [ self.matrix[i][j] for i in range(12) ]


class SetClass:
    '''
    Defines a set class in a musical system.
    '''

    def __init__(self, v, MS=MusicalSystem()):

        self.initv = list(v)
        self.rho = MS.rho
        self.phi = MS.phi
        self.MusicalSystem = MS

        self.normal_form = MS.normal_form(v)

        # now, the normal form is already the most compressed.
        # what we need it get the form that starts with the notes
        # "most to the right".

        w = self.normal_form
        self.prime_form = min((MS.vrho(MS.vphi(w, power=j), power=i)
                               for i in range(12) for j in range(2)),
                              key=lambda l: [MS.sveckey(x) for x in l])
        self.vec = list(self.prime_form)

    def __eq__(self, other):
        if not isinstance(other, SetClass):
            return False

        if self.MusicalSystem != other.MusicalSystem:
            return False

        return self.prime_form == other.prime_form

    def __repr__(self):
        return f"SetClass({self.initv}), with rho = {self.rho}, phi = {self.phi}.  Prime form: {print_prime_form(self.prime_form)}"

    def __str__(self):
        # return str(self.vec)
        return print_prime_form(self.prime_form)

    def int_vector(self):
        '''
        Returns the interval class vector with respect to rhoh.
        '''
        MS = self.MusicalSystem
        v = self.vec
        l = len(v)
        if MS.is_transitive:
            res = [0 for i in range(6)]
        else:
            res = [0 for i in range(9)]
        for i in range(l):
            for j in range(i+1, l):
                x = MS.interval_class(v[i], v[j])
                res[x-1] += 1
        return res

    def symmetries(self):
        '''
        Gives the powers of rho and phi that map v to w as sets.
        If only v is given, check inner symmetries.  (Always has
        identity as one.)
        '''
        sv = list(self.vec)  # aleready sorted
        MS = self.MusicalSystem
        return MS.symmetries(sv)

    def nsym(self):
        '''
        Number of symmetries that take v to itself.

        The result is a pair.  The first entry is how many
        powers of rhoh take v to v.

        The second is how many of rhoh^i*phi take v to v.
        '''
        sv = list(self.vec)  # aleready sorted
        MS = self.MusicalSystem
        sym = MS.symmetries(sv)
        sum1 = sum(1 for v in sym if v[1] == 0)
        sum2 = sum(1 for v in sym if v[1] == 1)
        return [sum1, sum2]

    def complement(self):
        '''
        Given a set class in prime form v, gives the set class of
        the complement of v.
        '''
        v = self.vec
        MS = self.MusicalSystem

        res = []
        for i in range(12):
            if i not in v:
                res.append(i)
        return SetClass(res, MS)

    def elements(self):
        '''
        Give all elements in the set class
        '''

        MS = self.MusicalSystem
        v = self.prime_form

        res = set()

        for i in range(12):
            for j in range(2):
                w = MS.vrho(MS.vphi(v, j), i)
                res.add(tuple(MS.normal_form(w)))

        return sorted(list(res))

    def class_sum(self):
        v = self.vec
        MS = self.MusicalSystem
        return MS.class_sum(v)

class ForteTable:
    '''
    Given n, gives the Forte Table of n-chords with respect
    to rho and phi.
    '''

    def __init__(self, n, MS=MusicalSystem()):

        # oreder with smaller set class size first
        if n > 6:
            n = 12 - 6

        self.n = n
        self.MusicalSystem = MS

        self.rho = MS.rho
        self.phi = MS.phi

        chords = set()
        if self.MusicalSystem.is_transitive:
            # a bit simpler: 0 is always first!
            for v in combinations(range(1, 12), n-1):
                w = [0] + list(v)
                s = tuple(SetClass(w, MS).vec)
                if not s in chords:
                    chords.add(s)
        else:
            for v in combinations(range(0, 12), n):
                w = list(v)
                s = tuple(SetClass(w, MS).vec)
                if not s in chords:
                    chords.add(s)
        chords = sorted(chords)

        res = []
        for v in chords:
            sc = SetClass(v, MS)
            scc = sc.complement()
            isc = sc.int_vector()
            iscc = scc.int_vector()
            s = sc.nsym()
            res.append([sc.vec, isc, s, iscc, scc.vec])

        self.matrix = res

    def __repr__(self):
        # return f"Forte Table for {self.n} and {12 - self.n} chords."
        return f"ForteTable({self.n})"

    def __str__(self):
        ft = self.matrix

        res = ''
        for v in ft:
            sc = print_prime_form(v[0])
            scc = print_prime_form(v[4])
            isc = vec_to_string(v[1], delim='[]')
            iscc = vec_to_string(v[3], delim='[]')
            sv = [sc, isc, '{:>2}, {:>2}'.format(v[2][0], v[2][1]), iscc, scc]
            res += '  '.join(sv) + '\n'

        res = res[:-1]

        return res

    def latex(self):
        '''
        Prints our a LaTeX formatted table.
        '''

        def chord_name(m):
            '''
            Name the chords.
            '''
            if m == 3:
                return 'Trichords'
            elif m == 4:
                return 'Tetrachords'
            elif m == 5:
                return 'Pentachords'
            elif m == 6:
                return 'Hexachords'
            elif m == 7:
                return 'Septachords'
            elif m == 8:
                return 'Octachords'
            elif m == 9:
                return 'Nonachords'

        ft = self.matrix
        n = self.n

        res = '\\begin{tabular}{llrll}\n'
        if n in [3, 4, 5, 6]:
            res += f'  \\multicolumn{{2}}{{c}}{{\\textbf{{{chord_name(n)}}}}} &  & \\multicolumn{{2}}{{c}}{{\\textbf{{{chord_name(12-n)}}}}} \\\\ \n'
            res += '  \\midrule \n'

        for v in ft:
            sc = '  $' + print_prime_form(v[0]) + '$'
            scc = '$' + print_prime_form(v[4]) + '$'
            # isc = '$' + vec_to_string(v[1], delim='[]') + '$'
            # iscc = '$' + vec_to_string(v[3], delim='[]') + '$'
            isc = '$' + ''.join([ str(x) for x in v[1] ]) + '$'
            iscc = '$' + ''.join([ str(x) for x in v[3] ]) + '$'
            sv = [sc, isc, f'${v[2][0]}$, ${v[2][1]}$', iscc, scc]
            res += ' & '.join(sv) + '\\\\ \n'

        res = res[:-4] + '\n'
        res += '\\end{tabular}'

        return res


#########################################################################

class ClassElements:
    '''
    Given n, it lists the elements in each class for n-chords.  At the "top"
    is the prime form, followed by the other elements in the same class.
    '''

    def __init__(self, n, MS=MusicalSystem()):

        self.rho = MS.rho
        self.phi = MS.phi
        self.n = n

        res = []
        if MS.is_transitive:
            for v in combinations(range(1, 12), n-1):
                vv = [0] + list(v)
                vv = tuple(SetClass(vv, MS).prime_form)
                for x in res:
                    if vv in x:
                        break
                else:
                    nv = SetClass(vv, MS).elements()
                    res.append(nv)
        else:
            for v in combinations(range(0, 12), n):
                vv = list(v)
                vv = SetClass(vv, MS).prime_form
                for x in res:
                    if vv in x:
                        break
                else:
                    nv = SetClass(vv, MS).elements()
                    res.append(nv)

        self.matrix = res

    def __repr__(self):
        return f'ClassElements({self.n})'

    def __str__(self):
        M = self.matrix
        n = self.n

        width = 2+2*(n-1)+2*n

        res = ''

        maxl = max([len(x) for x in M])

        for i in range(maxl):
            vstr = []
            for x in M:
                if i < len(x):
                    mstr = '{c:>{w}}'.format(c=vec_to_string(x[i]), w=width)
                else:
                    mstr = ' '*width
                vstr.append(mstr)
            res += '  '.join(vstr)
            res += '\n'
        res = res[:-1]

        return res

    def print_col(self, ncol=10):
        '''
        Prints the table of set classes breaking into ncol columns.
        '''
        M = self.matrix
        n = self.n

        width = 2+2*(n-1)+2*n

        l = len(M)
        res = ''

        q = l//ncol

        MM = [M[ncol*i: ncol*(i + 1)] for i in range(q)]
        if (l % ncol) != 0:
            MM += [M[ncol*q:]]

        for Mi in MM:
            maxl = max([len(x) for x in Mi])

            for i in range(maxl):
                res = ''
                vstr = []
                for x in Mi:
                    if i < len(x):
                        mstr = '{c:>{w}}'.format(
                            c=vec_to_string(x[i]), w=width)
                    else:
                        mstr = ' '*width
                    vstr.append(mstr)
                res += '  '.join(vstr)
                print(res)

            print('-------------------------------------------------------------\n')
