#!/usr/bin/env sage
# -*- mode: python -*-
"""
Hecke algebras on the homology of congruence subgroups of Euclidean Bianchi groups.

For Gamma_0(a) inside SL_2(O_d) / PSL / GL_2(O_d) / PGL, with d one of the five
Euclidean discriminants -1, -2, -3, -7, -11 and a = pi^e a prime power, this
decomposes H_1 into simultaneous eigenspaces for the Hecke operators and prints
the eigenvalue systems.

Coefficients are either
  * the free part, with eigenvalues in number fields -- exact arithmetic, no
    floating point, so an eigenvalue is reported as an element of an explicit
    field rather than a decimal; or
  * mod p for a rational prime p, with eigenvalues in GF(p^k).

Usage:
    sage bianchi_hecke.sage

The code needs Sage's preparser (it uses ^ for exponentiation), so run it with
`sage`, not with plain python.

Version:          2026-08-17
Generated from:   "New attempt at setup.ipynb", definition cells 0-3, copied verbatim
Core fingerprint: sha256:b3c6f2b3ff6e
"""

SCRIPT_VERSION   = '2026-08-17'
CORE_FINGERPRINT = 'b3c6f2b3ff6e'

import copy

def Norm(number):
    #alpha is a squarefree integer, defining the quadratic number field "Q adjoin sqrt(alpha)"
    #number is a list with two integers a and b, which determine an element inside the ring of integers of "Q adjoin sqrt(alpha)", where the first integer in the list is the integer component.
    #The function returns the absolute norm of the number.
    if alpha % 4 == 1:
        output = number[0]^2+number[0]*number[1]+ ((1-alpha)//4)*number[1]^2
    else:
        output = number[0]^2 - alpha * number[1]^2
    return output
def Is_Prime():
    #Checks whether prime_base is actually prime in the number field we're working in. Input is the prime_base 
    Num = list(factor(O(prime_base)))
    #print(Num)
    if len(Num) == 1 and Num[0][1] == 1:
        return True
    else:
        return False
def Split_Inert_Ramified():

    #Checks whether prime_base is a split, inert, or ramified prime over the ring of integers of the euclidian quadratic imaginary field
    #Outputs 0 if Split, 1 if Inert, and 2 if Ramified
    #The nature of the prime affects how right-coset representatives are constructed
    dummy = f(prime_base)
    Norm_dummy = Norm(dummy)
    #print(dummy, Norm_dummy)
    if dummy[1] == 0:
        return 1
    elif (alpha == -1 and Norm_dummy == 2) or (alpha == -Norm_dummy):
        return 2
    else:
        return 0
def Two_Power():
    #Checks whether level is a power of 2, returns True if it is, False if it isn't
    #Again, this affects how the Transversals are made
    Even_or_Odd = Norm(f(prime_base))
    if Even_or_Odd % 2 == 0:
        return True
    else:
        return False

def row_column_remover(matrix1):
    nraa = matrix1.nrows()
    nr2 = matrix1.ncols()
    #print(nraa, 'one')
    zeroes = vector([0] *nr2)
    #print(zeroes)# number of rows (this bit of code taken from https://ask.sagemath.org/question/47418/removing-the-repeated-rows-of-a-matrix/)
    eliminate = []# this will store the numbers of the rows to eliminate
    matrix2 = copy.deepcopy(matrix1)
    #matrix2 = matrix2.delete_columns([0])
    for i in range (nraa-1):# check all rows, except the last one
        #print(matrix2[i])
        if modulus(list(matrix2[i])) == 0:
            eliminate.append(i)
        for j in range(i+1, nraa):# compare to the rows below row i
            if matrix2[i] == matrix2[j] or matrix2[j] == -matrix2[i]:# or matrix1[i] == zeroes:
                if j not in eliminate:
                    eliminate.append(j)# if row j equals row i, list j to eliminate
    #print(matrix1.nrows(),nraa, 'twice')
    if nraa > 0:
        if matrix1[nraa-1] == zeroes:
            eliminate.append(nraa-1)
    matrix1 = matrix1.delete_rows(eliminate)
    return matrix1# eliminate the listed rows and return the resulting submatrix

def structure(smith_norm_matrix):
    mat = smith_norm_matrix.transpose()
    nr = mat.nrows()
    structure = [0]
    for i in range(nr):
        m = sum(mat[i])
        if m== 0:
            structure[0]+=1
        elif m > 1:
            structure.append(m)
    return structure

def modulus(numberlist):
    mod = 0
    for i in range(len(numberlist)):
        mod += abs(numberlist[i])^2
    return mod

def smallest(matrix1, colj):
    numberlist = matrix1.column(colj)
    start = True
    temp = 0
    m = 0
    for i in range(len(numberlist)):
        if numberlist[i] != 0 and start == True:
            temp = abs(numberlist[i])
            index = i
            start = False
            m = modulus(matrix1[i])
        elif abs(numberlist[i]) < temp and numberlist[i] != 0:
            temp = abs(numberlist[i])
            index = i
            m = modulus(matrix1[i])
        elif abs(numberlist[i]) == temp and modulus(matrix1[i]) < m:
            index = i
            m = modulus(matrix1[i])
    #print(numberlist[index], 'aaaa')
    return index

def fieldconstructor(alpha):
    x = polygen(ZZ, 'x')
    if alpha % 4 == 1:
        field.<x> = NumberField(x^2-x+(1-alpha)//4)
    else:
        field.<x> = NumberField(x^2 - alpha)
        #print(x^2)
    return x
    
def concatenate_words(word):
    newword = []
    negative_identity_counter = 0
    lastadder = True
    if word[0] == [0,0,2]:
        negative_identity_counter += 1
        del(word[0])
    i = len(word) - 1
    while i > 0:
        #print(i)
        while word[i][2] > 1:
            negative_identity_counter += 1
            word[i][2] -= 2
        if word[i] == [0,0,0]:
            i -= 1
            continue
        elif word[i] == [0,0,1]:
            new = [word[i-1][0],word[i-1][1],(word[i-1][2]+word[i][2]) % 4]
            if i == 1:
                lastadder = False
            i -= 1
        elif word[i-1][2] % 2 == 0:
            j = i - 1
            new = new = [word[i][0],word[i][1],word[i][2]]
            while word[j][2] % 2 == 0 and j > -1:
                new[0] += word[j][0]
                new[1] += word[j][1]
                new[2] += word[j][2]
                new[2] = new[2] % 4
                j-= 1
            i = j+1
            if i in [0, 1] and word[0][2] % 2 == 0:
                lastadder = False
        else:
            new = word[i]
        while new[2] > 1:
            negative_identity_counter +=1
            new[2] -= 2
        newword.insert(0, new)
        #print(newword)
        i -=1
    if lastadder == True and len(word) > 0:
        newword.insert(0, word[0])
    if negative_identity_counter % 2 == 1:
        if len(newword) > 0:
            newword[len(newword)-1][2] += 2
        else:
            newword = [[0,0,2]]
    return newword

def concatenate_loop(word):
    V = concatenate_words(word)
    W = concatenate_words(V)
    while len(V) > len(W):
        #print(W)
        temp = W
        W =concatenate_words(W)
        V = temp
        del(temp)
    return W

def word_to_vector(word):
    M = Matrix([[1,0],[0,1]])
    for i in range(len(word)):
        M = M*gen1^(word[i][0])*gen2^(word[i][1])*S^(word[i][2])
    return M

def Schreier_words(input_matrix_list):
    "Turns specific matrices that occur in the Schreier Transversals, into words in generators of Euclidian Bianchi groups"

    "This function does not work for arbitrary matrices"
    wordlist = []
    for i in range(len(input_matrix_list)):
        input_matrix = input_matrix_list[i]
        if input_matrix[1][1] == 1:
            z1 = input_matrix[1][0]
            word = [[f(z1)[0],f(z1)[1], 0]]
        elif input_matrix[1][0] == 1:
            z1 = -input_matrix[1][1]
            word = [[f(z1)[0],f(z1)[1], 1]]
        else:
            z1 = -input_matrix[1][1]
            z2 = -input_matrix[0][0]
            word = [[f(z1)[0],f(z1)[1],1],[f(z2)[0],f(z2)[1],0]]
        assert(word_to_vector(word) == input_matrix)
        wordlist.append(word)
    #print(wordlist, 'ha')
    return wordlist

def minimalquotient():#For making the  quotientring in general. This is very slow, and will only be used in the case of ramified prime base.
    ideal = O.fractional_ideal(level)
    quotientring = QuotientRing(O, level, names = 'a')
    number = O(level).absolute_norm()
    #print(number)
    lists = [1,-1, O(x), O(-1*x)]
    returnlist = [0]
    checklist = [quotientring(0)]
    condition = True
    while condition:
        for j in returnlist:
            for i in lists:
                #print(i,j)
                if (i+j not in returnlist) and (quotientring(i+j) not in checklist):
                    returnlist.append(i+j)
                    checklist.append(quotientring(i+j))
                if len(returnlist) >= number:
                    condition = False
                    break
    return returnlist
def Quotient_by_Level():
    "This function generates a quotient of the ring of integers O_K by the ideal generated by the level of the congruence subgroup. This procedure varies depending on the nature of the prime base of the level"
    sir = Split_Inert_Ramified()
    level_norm = Norm(f(level))
    prime_norm = Norm(f(prime_base))
    if sir == 0:
        even = Two_Power()
        if even == False:
            edge = (level_norm - 1) // 2
            returnlist = list(range(- edge, edge + 1))
        else:
            edge = level_norm // 2
            returnlist = list(range(- edge + 1, edge + 1))
    elif sir == 1:
        returnlist = []
        even = Two_Power()
        if even == False:
            edge = (level - 1) // 2
            candidates = list(range(- edge, edge + 1))
            for a in candidates:
                for b in candidates:
                    returnlist.append(a+b*x)
        else:
            #Residues here are a+b*x with a and b each running over a complete
            #residue system mod level, so there must be `level` candidates, not
            #level_norm.  This matches the odd branch above, which uses
            #(level-1)//2.  Only reachable for an inert prime 2.
            edge = level // 2
            candidates = list(range(- edge + 1, edge + 1))
            for a in candidates:
                for b in candidates:
                    returnlist.append(a+b*x)
    else:
        returnlist = minimalquotient()
    return returnlist
def ZeroDivsQuotient():
    zero_divisors = []
    for number in Fixed_Quotient:
        if number != 0 and Norm(f(number)) != 1:
            #print(number,  gcd(O(level),O(number)))
            if gcd(O(level), O(number)) != 1:
                zero_divisors.append(number)
    return zero_divisors
def Schreier_Transversal_New():
    level_norm = Norm(f(level))
    prime_norm = Norm(f(prime_base))
    Acosets = []
    Scosets = []
    Input_list = Quotient_by_Level()
    if Split_Inert_Ramified() == 1:
        base_norm = prime_base
    else:
        base_norm = prime_base.absolute_norm()
    for item in Input_list:
        M = matrix([[1,0],[item, 1]])
        Acosets.append(M)
        if O(item).absolute_norm() % base_norm == 0:
            Scosets.append(M*S)
    cosets = Acosets + Scosets
    Checker = level_norm + level_norm // prime_norm
    #print(Checker, len(cosets))
    assert(len(cosets) == Checker)
    return cosets  
    
def invert_word(word):
    if word == [[0,0,0]] or word == [[0,0,2]]:
        return word
    newword = []
    for i in reversed(range(0,len(word))):
        if word[i][0] == word[i][1] == 0:
            new = [[0,0,-word[i][2] %4]]
        elif word[i][2] != 0:
            new = [[0,0,-word[i][2] % 4],[-word[i][0],-word[i][1],0]]
        else:
            new = [[-word[i][0],-word[i][1],0]]
        newword = newword +new
    #print(newword)
    return concatenate_loop(newword)

def power_word(word, exponent):
    if exponent == 0:
        return [[0,0,0]]
    newword = word.copy()
    for i in range(1,exponent):
        newword += word
    newword = concatenate_loop(newword)
    return newword

def Take_out_negative_identity(gens):
    #print(gens, 'here')
    if gens != [[0,0,2]]:
        for bit in gens:
            bit[2] = bit[2] % 2
    #print(gens, 'there')
    return gens

def Subgroup_generators():
    full_group_generators = [gen1, gen2, S]
    #print('take long?',O)
    Cosets_Matrices = Schreier_Transversal_New()
    Cosets_words = Schreier_words(Cosets_Matrices)
    #A coset is a point [a:b] of P^1(O/level).  The level is always a prime
    #power, and a determinant-1 matrix has (a,b) generating the unit ideal, so
    #the prime cannot divide both: at least one of a,b is invertible mod level.
    #That gives every point a unique normal form -- [a/b : 1] when b is a unit,
    #otherwise [1 : b/a] -- which is unchanged if (a,b) is rescaled by a unit.
    #I.reduce() is a canonical representative of a residue class, so the normal
    #form becomes a hashable key and the matching coset can be looked up in a
    #dict instead of being compared against every representative in turn.
    I = O.fractional_ideal(level)
    unit_ideal = O.fractional_ideal(O(1))
    def unit_mod_level(z):
        z = O(z)
        if z == 0:
            return False
        return I + O.fractional_ideal(z) == unit_ideal
    def coset_key(matr):
        a = O(matr[1][0]); b = O(matr[1][1])
        if unit_mod_level(b):
            return ('A', tuple(f(O(I.reduce(a*b.inverse_mod(I))))))
        assert unit_mod_level(a), 'coset point is not primitive - is the level a prime power?'
        return ('S', tuple(f(O(I.reduce(b*a.inverse_mod(I))))))
    Coset_index = {}
    for j in range(len(Cosets_Matrices)):
        Coset_index[coset_key(Cosets_Matrices[j])] = j
    Returnlist = []
    for mat in full_group_generators:
        for i in range(len(Cosets_Matrices)):
            j = Coset_index[coset_key(Cosets_Matrices[i]*mat)]
            Gen_word = concatenate_loop(Cosets_words[i] + Schreier_words([mat])[0]+ invert_word(Cosets_words[j]))
            Gen_word2 = concatenate_loop(invert_word(Gen_word)+[[0,0,2]])
            if Gen_word != [[0,0,0]] and  Gen_word not in Returnlist and invert_word(Gen_word) not in Returnlist and Gen_word2 not in Returnlist:
                temp = word_to_vector(Gen_word)[1][0]
                assert((temp*(level^(-1))).denominator() == 1)
                Returnlist.append(Take_out_negative_identity(Gen_word))
    return Returnlist

def special_generators():
    Gens = Subgroup_generators()
    #print(Gens, "WAWOWIWA")
    number_ring = O
    #p =x.coordinates_in_terms_of_powers()
    Gens.remove([[0,0,2]])
    special_gens = [[[0,0,2]]]
    remover = []
    L = f(O(level))
    special_gens.append([[L[0],L[1],0]])
    L2 = (f(O(level*x)))
    special_gens.append([[L2[0],L2[1],0]])

    for i in Gens:
        if len(i) == 1:
            remover.append(i)
    for j in Gens:
        if len(j) == 2 and j[0] == [0,0,1]:
            remover.append(j)
    for y in remover:
        Gens.remove(y)
    special_gens.append([[0,0,1],[1,0,3]])
    special_gens.append([[0,0,1],[0,1,3]])
    return special_gens + Gens

def Matrix_Generators():
    TheList = []
    if alpha in [-1,-3]:
        TheList.append(Matrix([[x,0],[0,x^(-1)]]))
    else:
        TheList.append(S^2)
    for index in range(1, len(Real_Gens)):
        TheList.append(word_to_vector(Real_Gens[index]))
    return TheList
    
Modrel = [[[-1,0,1],[-1,0,1],[-1,0,3]]]
Bianchirel1 = [[[0,-1,1],[0,1,1],[0,-1,1],[0,1,1],[0,-1,1],[0,1,3]],
              [[0,-2,1],[0,1,1],[0,-2,1],[0,1,3]],
              [[-1,-1,1],[-1,1,1],[-1,-1,1],[-1,1,3]],
              ]
Bianchirel2 = [[[0,1,1],[0,-1,1],[0,1,1],[0,-1,3]]]
Bianchirel3 = [[[0,0,3],[0,-1,1],[-2,2,1],[0,-1,1],[-2,2,0]],
                [[0,0,1],[-1,-1,1],[-1,1,1],[1,-1,1],[1,-1,1],[-1,1,0]],
              [[0,0,1],[0,-1,1],[-1,1,1],[0,-1,1],[-1,1,1],[0,-1,1],[-1,1,0]]]
Bianchirel7 = [[[-1,1,1],[0,-1,1],[-1,1,1],[0,-1,3]]]
Bianchirel11 = [[[-1,1,3],[0,-1,1],[-1,1,1],[0,-1,1],[-1,1,1],[0,-1,1]]]

def get_relations():
    if alpha == -1:
        return copy.deepcopy(Modrel +Bianchirel1)
    if alpha == -2:
        return copy.deepcopy(Modrel +Bianchirel2)
    if alpha == -3:
        return copy.deepcopy(Modrel +Bianchirel3)
    if alpha == -7:
        return copy.deepcopy(Modrel +Bianchirel7)
    if alpha == -11:
        return copy.deepcopy(Modrel +Bianchirel11)
    else:
        sys.exit()


def Unit_Relation():
    if alpha == -1:
        returnrel = [[0,1,1],[0,-1,1],[0,1,1]]
    elif alpha == -3:
        returnrel = [[0,1,1],[1,-1,1],[0,1,1]]
    else:
        returnrel = [[-1,0,1],[-1,0,1],[-1,0,1]]

    FinalRel = Rewrite_Relation(returnrel)
    return FinalRel
    
def Reidemeister_Schreier_Relations():

    #relations2 = copy.deepcopy(relations)
    relations = get_relations()
    Transversal = Schreier_words(Schreier_Transversal_New())
    Returnlist = []
    for j in Transversal:
        if j == [[0,0,0]]:
            for i in relations:
                Returnlist.append(copy.deepcopy(i))
        else:
            for i in relations:
                #print(j)
                #print(j, invert_word(j))
                newword = concatenate_loop(j + i + invert_word(j))
                Returnlist.append(newword)
    #print(Returnlist, 'ziziziz')
    return Returnlist


def Quotient_Match(relationstring,relationword):
    number = relationstring[0][0]+relationstring[0][1]*x
    for item in Fixed_Quotient:
        divided = (number-item)*(level^(-1))
        if divided.denominator() == 1:
            num = f(divided)
            relationword.append([1,num[0]])
            relationword.append([2,num[1]])
            relationstring[0][0] = f(item)[0]
            relationstring[0][1] = f(item)[1]
            break
    return [relationstring, relationword]
def Seperate_Inv_generators(generators):
    lengthgens = len(generators)
    rels = []
    for i in range(lengthgens):
        if i > 4 and len(generators[i]) == 2:
            rels.append([generators[i],i])
    return rels
def Seperate_PriPow_generators(generators):
    lengthgens = len(generators)
    real_rels = []
    im_rels = []
    for i in range(lengthgens):
        if len(generators[i]) > 2:
            if generators[i][1][0] == 1:
                real_rels.append([generators[i],i])
            else:
                im_rels.append([generators[i],i])
    if len(real_rels) > 0 and len(im_rels) > 0:
        return [real_rels, im_rels]
    else:
        return []

def PriPow_relation(Genies, relate_string, relation):


    #THIS ONLY WORKS FOR PRIME-POWER LEVELS
    #Input: "generators" is the output of "Seperate_PriPow_generators"
    #Input: "vector" is the vector associated to the relation generated by Reidemeister-Schreier
    #Input: "Relation" is the string that is being written in terms of generators for the relevant congruence subgroup
    #print('This is being read', relation)
    endrelator = []
    real_gens = Genies[0]
    real_len = len(real_gens)
    im_gens = Genies[1]
    im_len = len(im_gens)
    num1 = relation[1][0]
    temp1 = 0
    num2 = relation[1][1]
    temp2 = 0
    startrel = relation[0]
    if num1 > 0:
        while temp1 != num1:
            for i in range(real_len):
                if real_gens[i][0][0] == startrel:
                    temp1 += 1
                    startrel = [-real_gens[i][0][2][0], -real_gens[i][0][2][1], 1]
                    relate_string.append([real_gens[i][1],1])
                    #shared with Real_Gens; concatenate_loop below mutates it
                    endrelator +=copy.deepcopy(real_gens[i][0])
                    break
    elif num1 < 0:
         while temp1 != num1:
            for i in range(real_len):
                if  startrel == [-real_gens[i][0][2][0], -real_gens[i][0][2][1], 1]:
                    temp1 -= 1
                    startrel = real_gens[i][0][0]
                    relate_string.append([real_gens[i][1],-1])
                    endrelator+=invert_word(real_gens[i][0])
                    break
    if num2 > 0:
        while temp2 != num2:
            for i in range(im_len):
                if im_gens[i][0][0] == startrel:
                    temp2 += 1
                    startrel = [-im_gens[i][0][2][0], -im_gens[i][0][2][1], 1]
                    relate_string.append([im_gens[i][1],1])
                    endrelator+=copy.deepcopy(im_gens[i][0])
                    break
    elif num2 < 0:
         while temp2 != num2:
            for i in range(im_len):
                if  startrel == [-im_gens[i][0][2][0], -im_gens[i][0][2][1], 1]:
                    temp2 -= 1
                    startrel = im_gens[i][0][0]
                    relate_string.append([im_gens[i][1],-1])
                    endrelator+=invert_word(im_gens[i][0])
                    break
    #print(endrelator)
    endrelator = concatenate_loop(endrelator)
    #print('This is what it made', endrelator)
    #print(endrelator[0], relation[0], endrelator[1], relation[1])
    assert(endrelator[0] == relation[0] and endrelator[1] == relation[1])
    relation = concatenate_loop(invert_word(endrelator) + relation)
    return [relation, relate_string]

def norm_list(listi):
    returnlist = []
    for i in range(5,len(listi)):
        if [listi[i][0][0],listi[i][0][1]] not in returnlist:
            returnlist.append([listi[i][0][0],listi[i][0][1]])
        if [-listi[i][len(listi[i])-1][0],-listi[i][len(listi[i])-1][1]] not in returnlist:
            returnlist.append([-listi[i][len(listi[i])-1][0],-listi[i][len(listi[i])-1][1]])
    return returnlist

def Rewriting_process():
    if alpha == -1:
        half_unit_ord = 2
    elif alpha == -3:
        half_unit_ord = 3
    else:
        half_unit_ord = 1
    I = O.fractional_ideal(level)
    generators = Real_Gens
    inv_gens = []
    zerodiv_gens = Seperate_PriPow_generators(Real_Gens)
    #print('Hoo Ha', zerodiv_gens)
    for l in range(len(generators)):
        #print(l, generators[l])
        if l > 4 and len(generators[l]) == 2:
            inv_gens.append([generators[l],l])
    relations =  Real_Rels
    quotientlist = Fixed_Quotient
    zero_divisors = []
    for number in quotientlist:
        if number != 0 and Norm(f(number)) != 1:
            if (level*number^(-1)).denominator() == 1:
                zero_divisors.append(number)
    #print(zero_divisors)
    norm = norm_list(generators)
    #print(quotientlist)
    #print(norm, 'hello')
    #print(number_ring)
    new_relations = []
    for gen in inv_gens:
        new_rel = []
        if gen[0][0][0] + gen[0][0][1] * x + gen[0][1][0] + gen[0][1][1] * x == 0:
            new_rel.append([0,half_unit_ord])
            new_rel.append([gen[1],2])
            new_relations.append(new_rel)
        elif gen[0][0][0] + gen[0][0][1] * x + gen[0][1][0] + gen[0][1][1] * x == 1:
            new_rel.append([gen[1],3])
            new_relations.append(new_rel)
        elif gen[0][0][0] + gen[0][0][1] * x + gen[0][1][0] + gen[0][1][1] * x == -1:
            new_rel.append([gen[1],3])
            new_rel.append([0,half_unit_ord])
            new_relations.append(new_rel)
    for rel in relations:
        new_rel = []
        word = copy.deepcopy(rel)
        condition = True
        while condition:
            #print(word, "AAA")
            word[0][2] = word[0][2] % 4
            if word[0] == [0,0,2]:
                word = concatenate_loop(generators[0]+word)
                new_rel.insert(0, [0,half_unit_ord])
            elif word[0] == [0,0,1]:
                if not word[1][0] == 0:
                    new_rel.append([3, word[1][0]])
                if not word[1][1] == 0:
                    new_rel.append([4, word[1][1]])
                word = concatenate_loop([[0,0,1],[-word[1][0],-word[1][1],3]]+word)
            elif word[0][0]+word[0][1]*x not in Fixed_Quotient and word[0] != [0,0,0]:
                Thing = Quotient_Match(word,new_rel)
                word = Thing[0]
                new_rel = Thing[1]
                #z1 = word[0][0]+word[0][1]*x
                #for i in quotientlist:
                    #if z1 - i in I:
                        #print(z1, i, z1-i)
                        #z2 = (z1-i)*level^(-1)
                        #L = f(z2)
                        #print(z2, L)
                        #break
                #if L[0] > 0:
                #    word = concatenate_loop(invert_word(power_word(generators[1], L[0])) + word)
                #elif L[0] < 0:
                #    word = concatenate_loop(power_word(generators[1],-L[0])+word)
                #if L[1] > 0:
                #    word = concatenate_loop(invert_word(power_word(generators[2], L[1]))+ word)
               # elif L[1] < 0:
                #    word = concatenate_loop(power_word(generators[2],-L[1])+word)
               # if L[0] != 0:
               #     new_rel.append([1, L[0]])
               # if L[1] != 0:
                #    new_rel.append([2,L[1]])
            elif word[0][2] % 2 == 1 and  word[0][0] + word[0][1] * x not in zero_divisors:
                for gen in inv_gens:
                    if gen[0][0][0] == word[0][0] and gen[0][0][1] == word[0][1]:
                        new_rel.append([gen[1],1])
                        word = concatenate_loop(invert_word(gen[0])+word)
                    elif gen[0][1][0] == -word[0][0] and gen[0][1][1] == -word[0][1]:
                        new_rel.append([gen[1],-1])
                        word = concatenate_loop(gen[0] + word)
            elif word[0][2] % 2 == 1 and  word[0][0] + word[0][1] * x in zero_divisors:
                temp = PriPow_relation(zerodiv_gens, new_rel, word)
                word = temp[0]
                new_rel = temp[1]
            else:
                word[0][2] = word[0][2] % 4
                assert(word == [[0,0,0]])
                condition = False
        #print(new_rel)
        M = Matrix([[1,0],[0,1]])
        for l in range(len(new_rel)):
            M = M*(word_to_vector(generators[new_rel[l][0]])^(new_rel[l][1]))
        #print(M, new_rel)
        if M == Matrix([[-1,0],[0,-1]]):
            if new_rel[0] == [0,half_unit_ord]:
                new_rel.remove([0,half_unit_ord])
            else:
                new_rel.insert(0,[0,half_unit_ord])
        else:
            assert(M == Matrix([[1,0],[0,1]]))
        new_relations.append(new_rel)
        #print(new_relations)
    return [new_relations, generators]


def SAS_relator(relation, new_relator):

    new_rel = new_relator
    word =relation 
    if not word[1][0] == 0:
        new_rel.append([3, word[1][0]])
    if not word[1][1] == 0:
        new_rel.append([4, word[1][1]])
    word = concatenate_loop([[0,0,1],[-word[1][0],-word[1][1],3]]+word)
    return [word, new_rel]

def Fin_Oder_Inv_gen():
    new_rels = []
    for gen in Inv_gens:
        num1 = gen[0][0][0] + gen[0][0][1]*x + gen[0][1][0] + gen[0][1][1]*x
        #num1 is minus the trace, so the characteristic polynomial is
        #lambda^2 + num1*lambda + 1 and the generator's order follows from it:
        #  num1 =  0  ->  lambda^2 + 1      -> order 4, gen^2 = -I
        #  num1 =  1  ->  lambda^2+lambda+1 -> order 3, gen^3 =  I
        #  num1 = -1  ->  lambda^2-lambda+1 -> order 6, gen^3 = -I
        #Only the order-4 and order-6 cases pick up the -I, recorded as
        #generator 0 to the power half_unit_ord.
        if num1 == 0:
            new_rels.append([[0,half_unit_ord ],[gen[1],2]])
        elif num1 == 1:
            new_rels.append([[gen[1],3]])
        elif num1 == -1:
            new_rels.append([[0,half_unit_ord],[gen[1],3]])
    return new_rels
        
def GL_Relator_Strings():
    "This adds relations that are given by conjugation by the diagonal matrix (1,u) by the semidirect product construction of GL from SL, built in combinatorial group computations on how it affects chosen generators"
    The_Relators = []
    #print(Real_Gens)
    if alpha == -1:
        the_unit = x
        unit_rel = [[0,-1,1],[0,1,1],[0,-1,1]]
    elif alpha == -3:
        the_unit = x
        unit_rel = [[1,-1,1],[0,1,1],[1,-1,1]]
    else: 
        the_unit = -1
        unit_rel = [[0,0,2]]
    for index in [1,2]:
        #print(Real_Gens[index],"yayayay")
        temp = Real_Gens[index][0][0] + Real_Gens[index][0][1]*x
        #print(temp, the_unit*temp, "aaaaa")
        number = temp - temp*the_unit
        #print(number, "booboo")
        numvec = f(number)
        The_Relators.append([[numvec[0],numvec[1],0]])
    numvec1 = f(O(1)-O(the_unit^(-1)))
    numvec2 = f(O(x)-O(x*the_unit^(-1)))
    The_Relators.append([[0,0,1],[numvec1[0],numvec1[1],3]])
    The_Relators.append([[0,0,1],[numvec2[0],numvec2[1],3]])
    for genie in Inv_gens:
        ngenie = genie[0]
        #print(ngenie[0],ngenie[0][0], (ngenie[0][0] + ngenie[0][1]*x),(-the_unit^(-1)), "OOOWIIEE")
        num1 = f((ngenie[0][0] + ngenie[0][1]*x)*(-the_unit^(-1)))
        num2 = f((ngenie[1][0] +ngenie[1][1]*x)*(-the_unit))
        #concatenate_loop mutates its argument in place, so unit_rel must be
        #copied -- otherwise the first iteration rewrites [[0,0,2]] to [[0,0,0]]
        #and every later relator is built with an empty unit_rel.
        newrel = concatenate_loop(ngenie + [[num2[0],num2[1],3],[num1[0],num1[1],0]] + copy.deepcopy(unit_rel))
        The_Relators.append(newrel)
        #print(ZeroDiv_gens)
    if len(ZeroDiv_gens) > 0:
        #Seperate_PriPow_generators returns [real_rels, im_rels], a PAIR of
        #lists, so iterating over ZeroDiv_gens directly walks the two sublists
        #and genie[0][0] then picks only the first generator of each -- two
        #relators regardless of how many zero divisors there are.  Concatenate
        #the sublists so every generator contributes, as PriPow_relation does.
        for genie in ZeroDiv_gens[0] + ZeroDiv_gens[1]:
            ngenie = genie[0]
            num1 = f((ngenie[0][0] + ngenie[0][1]*x)*(-the_unit))
            num2 = f((ngenie[1][0] +ngenie[1][1]*x)*(-the_unit^(-1)))
            num3 = f((ngenie[2][0] + ngenie[2][1]*x)*(-the_unit))
            #Checking D g D^-1 against the inverse of this tail shows the two
            #differ by exactly -I for every zero-divisor generator, in every
            #field.  -I is central and diagonal matrices lie in every Gamma_0,
            #so the extra [[0,0,2]] can sit anywhere in the word.
            newrel = concatenate_loop(ngenie + [[num3[0],num3[1],1],[num2[0],num2[1],3],[num1[0],num1[1],0]] + [[0,0,2]])
            The_Relators.append(newrel)
    return The_Relators
def Rewrite_Relation(relation):
    new_rel = []
    if alpha == -1:
        half_unit_ord = 2
    elif alpha == -3:
        half_unit_ord = 3
    else:
        half_unit_ord = 1
    word = copy.deepcopy(relation)
    condition = True
    while condition:
        #print(word, "AAA")
        #print(Inv_gens, "YOOYOOYOO")
        word[0][2] = word[0][2] % 4
        #S^3 = -I * S, and -I is central, so pull the -I out into the unit
        #generator and leave S.  Without this a word starting [a,b,3] falls
        #through to the Inv_gens branch, which cannot rewrite it, and the loop
        #spins forever (hit by Bianchirel3[0] = [[0,0,3],...] for alpha = -3).
        if word[0][2] % 4 == 3:
            new_rel.insert(0, [0,half_unit_ord])
            word[0][2] = 1
        if word[0] == [0,0,0] and len(word) == 1:
            condition = True
            break
        elif word[0] == [0,0,2]:
            #concatenate_loop mutates its argument, and Real_Gens[0] is shared
            #state, so it must be copied or the generator list gets corrupted.
            word = concatenate_loop(copy.deepcopy(Real_Gens[0])+word)
            new_rel.insert(0, [0,half_unit_ord])
        elif word[0] == [0,0,1]:
            Thing = SAS_relator(word,new_rel)
            word = Thing[0]
            new_rel = Thing[1]
        elif word[0][0]+word[0][1]*x not in Fixed_Quotient:
            Thing = Quotient_Match(word,new_rel)
            word = Thing[0]
            new_rel = Thing[1]
        elif word[0][2] % 2 == 1 and  word[0][0] + word[0][1] * x not in zero_div:
            for gen in Inv_gens:
                #print(gen, "HALLOOO")
                if gen[0][0][0] == word[0][0] and gen[0][0][1] == word[0][1]:
                    new_rel.append([gen[1],1])
                    word = concatenate_loop(invert_word(gen[0])+word)
                elif gen[0][1][0] == -word[0][0] and gen[0][1][1] == -word[0][1]:
                    new_rel.append([gen[1],-1])
                    word = concatenate_loop(copy.deepcopy(gen[0]) + word)
        elif word[0][2] % 2 == 1 and  word[0][0] + word[0][1] * x in zero_div:
            temp = PriPow_relation(ZeroDiv_gens, new_rel, word)
            word = temp[0]
            new_rel = temp[1]
        else:
            word[0][2] = word[0][2] % 4
            assert(word == [[0,0,0]])
            condition = False
    #print(new_rel)
    #Generator 0 is the unit matrix diag(x, x^-1), of order 2*half_unit_ord with
    #gen0^half_unit_ord = -I.  word_to_vector(Real_Gens[0]) is -I, NOT that
    #matrix, so the check must be made against Real_Mat_Gens; otherwise a
    #relation that really equals -I gets recorded as trivial, which silently
    #imposes the PSL relation and halves the 2-part of the torsion.
    M = Matrix([[1,0],[0,1]])
    for l in range(len(new_rel)):
        M = M*(Real_Mat_Gens[new_rel[l][0]]^(new_rel[l][1]))
    unit_gen = Real_Mat_Gens[0]
    corrected = False
    for i in range(2*half_unit_ord):
        if M == unit_gen^i:
            if i != 0:
                new_rel.insert(0, [0, -i])
            corrected = True
            break
    assert corrected, 'rewritten relation is not a power of the unit generator'
    return new_rel
    
def Rewrite_Relation_GL(relation):
    new_rel = []
    #Unlike Rewrite_Relation, the input here is NOT a relator: it is an
    #element whose class is killed in the GL abelianisation.  So the test
    #is not 'does this rewrite to the identity' but 'does it reproduce the
    #original matrix'.  Record it now, before the rewriting consumes it.
    Original_matrix = word_to_vector(copy.deepcopy(relation))
    if alpha == -1:
        half_unit_ord = 2
    elif alpha == -3:
        half_unit_ord = 3
    else:
        half_unit_ord = 1
    word = copy.deepcopy(relation)
    condition = True
    while condition:
        #print(word, "AAA")
        #print(Inv_gens, "YOOYOOYOO")
        word[0][2] = word[0][2] % 4
        #S^3 = -I * S, and -I is central, so pull the -I out into the unit
        #generator and leave S.  Without this a word starting [a,b,3] falls
        #through to the Inv_gens branch, which cannot rewrite it, and the loop
        #spins forever (hit by Bianchirel3[0] = [[0,0,3],...] for alpha = -3).
        if word[0][2] % 4 == 3:
            new_rel.insert(0, [0,half_unit_ord])
            word[0][2] = 1
        if word[0] == [0,0,0] and len(word) == 1:
            condition = True
            break
        elif word[0] == [0,0,2]:
            #concatenate_loop mutates its argument, and Real_Gens[0] is shared
            #state, so it must be copied or the generator list gets corrupted.
            word = concatenate_loop(copy.deepcopy(Real_Gens[0])+word)
            new_rel.insert(0, [0,half_unit_ord])
        elif word[0] == [0,0,1]:
            Thing = SAS_relator(word,new_rel)
            word = Thing[0]
            new_rel = Thing[1]
        elif word[0][0]+word[0][1]*x not in Fixed_Quotient:
            Thing = Quotient_Match(word,new_rel)
            word = Thing[0]
            new_rel = Thing[1]
        elif word[0][2] % 2 == 1 and  word[0][0] + word[0][1] * x not in zero_div:
            for gen in Inv_gens:
                #print(gen, "HALLOOO")
                if gen[0][0][0] == word[0][0] and gen[0][0][1] == word[0][1]:
                    new_rel.append([gen[1],1])
                    word = concatenate_loop(invert_word(gen[0])+word)
                elif gen[0][1][0] == -word[0][0] and gen[0][1][1] == -word[0][1]:
                    new_rel.append([gen[1],-1])
                    word = concatenate_loop(copy.deepcopy(gen[0]) + word)
        elif word[0][2] % 2 == 1 and  word[0][0] + word[0][1] * x in zero_div:
            temp = PriPow_relation(ZeroDiv_gens, new_rel, word)
            word = temp[0]
            new_rel = temp[1]
        else:
            word[0][2] = word[0][2] % 4
            assert(word == [[0,0,0]])
            condition = False
    #print(new_rel)
    M = Matrix([[1,0],[0,1]])
    for l in range(len(new_rel)):
        M = M*(Real_Mat_Gens[new_rel[l][0]]^(new_rel[l][1]))
    #If the rewriting came out as gen0^i times the original instead of the
    #original itself, record the compensating gen0^-i.  gen0 is the diagonal
    #unit matrix, which lies in every Gamma_0, so this is pure bookkeeping.
    unit_gen = Real_Mat_Gens[0]
    for i in range(2*half_unit_ord):
        if M == (unit_gen^i)*Original_matrix:
            if i != 0:
                new_rel.insert(0, [0, -i])
            break
    return new_rel

def proj_equality(proj_space, vector1, vector2):
    x = proj_space([vector1[0],vector1[1]])
    y = proj_space([vector2[0],vector2[1]])
    if x == y:
        return True
    else:
        return False
def projpoint(proj_space, matr1):
    vector1 = vector(matr1[1,:])
    x = proj_space([vector1[0],vector1[1]])
    return x

def unit_matrix(unit):
    thismatrix= matrix(([unit,0],[0,unit^(-1)]))
    assert(thismatrix.determinant() == 1)
    return thismatrix

def alpha_coeff():
    if alpha % 4 == 1:
        return (1-alpha)//4
    else:
        return -alpha
def norms(number):
    if number in ZZ:
        return number^2
    else:
        return number.absolute_norm()

def Euclidian_Quotient(Higher_number, Lower_number):
    Approx_quotient = Higher_number*(Lower_number^(-1))
    Approx_Vec = f(Approx_quotient)
    real_part = [Approx_Vec[0].floor(),Approx_Vec[0].ceil()]
    im_part = [Approx_Vec[1].floor(),Approx_Vec[1].ceil()]
    #print(real_part,im_part)
    qu1 = O(real_part[0]) + O(im_part[0])*O(x)
    qu2 = O(real_part[1]) + O(im_part[0])*O(x)
    qu3 = O(real_part[0]) + O(im_part[1])*O(x)
    qu4 = O(real_part[1]) + O(im_part[1])*O(x)
    op1 = norms(Higher_number - qu1*Lower_number)
    op2 = norms(Higher_number - qu2*Lower_number)
    op3 = norms(Higher_number - qu3*Lower_number)
    op4 = norms(Higher_number - qu4*Lower_number)
    qulist = [qu1, qu2, qu3, qu4]
    oplist = [op1, op2, op3,op4]
    #print(qulist, "ya", oplist)
    min_index = oplist.index(min(oplist))
    return qulist[min_index]
    
def bounded_norm(number):#, field):
    if alpha >= 0:
        print('the function does not support that')
        sys.exit()
    bounded = []
    x = 0
    while ((x+1)^2 <= number):
        x += 1
    bounded += [x]
    y = 0
    while(alpha_coeff()*(y+1)^2 <= number):
        y +=1
    bounded += [y]
    circle = []
    highest_norm = 0
    options = []
    for z1 in range(-bounded[0], bounded[0]+1):
        for z2 in range(-bounded[1], bounded[1]+1):
            if alpha % 4 == 1:
                if z1^2 +z1*z2 +((1 - alpha)//4)*(z2^2) <= number:
                    circle += [[z1, z2]]
                    if z1^2 +z1*z2 +((1 - alpha)//4)*(z2^2) > highest_norm:
                        options = [[z1,z2]]
                        highest_norm  = z1^2 +z1*z2 +((1 - alpha)//4)*(z2^2)
                    elif z1^2 +z1*z2 +((1 - alpha)//4)*(z2^2) == highest_norm:
                        options += [[z1,z2]]
            else:
                if z1^2 -alpha*(z2^2) <= number:
                    circle += [[z1,z2]]
                    if z1^2 -alpha*(z2^2) > highest_norm:
                        options = [[z1,z2]]
                        highest_norm  = z1^2 -alpha*(z2^2)
                    elif z1^2 -alpha*(z2^2) == highest_norm:
                        options += [[z1,z2]]

    return circle

def Afunc(el_of_number_field):
    assert(denominator(el_of_number_field) == 1)
    Z = matrix( ([1,0],[el_of_number_field,1]))
    return Z

def remainder(cq_integer1,cq_integer2):
    norm1 = norms(cq_integer1)
    #print(norm1)
    norm2 = norms(cq_integer2)
    #print(norm2)
    options = bounded_norm(norm1)
    #print(options,"yayaya")
    thing = norm2
    t = 0
    for i in range(len(options)):
        r = -cq_integer2*(options[i][0]+options[i][1]*x)+cq_integer1
        if norms(r) < thing:
            output = r
            thing = norms(r)
            #print('ho',norms(r))
            #print('hi',r)
            t = -options[i][0]-options[i][1]*x
    return t

def Euclidian_algorithm(matric):

    full_relatorlist = []
    identitymatrix = Matrix(([1,0],[0,1]))
    Det = matric.determinant()
    #Determinant 1 means no GL twist was needed; every element of Gamma_0(level)
    #lands here, so GLFactor must be defined before the branch rather than only
    #inside it.
    GLFactor = 0
    if Det != 1:
        if alpha == -1:
            for i in range(1,4):
                if Det == x^i:
                    M = matric([[1,0],[0, Det^(-1)]])
                    matric = matric * M
            GLFactor = i
        elif alpha == -3:
            for i in range(1,6):
                if Det == x^i:
                    M = matric([[1,0],[0, Det^(-1)]])
                    matric = matric * M            
            GLFactor = i
        else:
            M = matric([[1,0],[0, Det]])
            matric = matric * M  
            GLFactor = 1 
    assert(matric.determinant() == 1)
    #Keep the (possibly GL-untwisted) input: the word must rebuild it exactly,
    #not merely up to the central -I.
    Target_matrix = copy.deepcopy(matric)
    relatorlist = [[0,0,0]]
    while(matric != identitymatrix ):
        #print(matric)
        if matric[0][1] == 0:
            if matric[0][0] == 1:
                L = f(matric[1][0])
                relatorlist.insert(0,[L[0],L[1],0])
                matric = matric*(Afunc(-matric[1][0]))
                continue
            elif matric[0][0] == -1:
                relatorlist[0][2] += 2
                matric = matric*unit_matrix(-1)
                assert(matric.determinant() == 1)
                continue
            else:
                unit = matric[0][0]
                unit2 = unit^(-1)
                u1 = f(unit)
                u2 = f(unit2)
                relator = [[ -u2[0],-u2[1],1],[-u1[0],-u1[1],1],[-u2[0],-u2[1],3]]
                relatorlist = relator +relatorlist
                matric = matric*unit_matrix(unit2)
                continue
        elif matric[0][0] == 0:
            matric = matric*S
            relatorlist.insert(0,[0,0,3])
            continue
        elif norms(matric[0][0]) == norms(matric[0][1]) and denominator(matric[0][0]*(matric[0][1])^(-1)) == 1:
            #Equal norms only give a single exact division step when the ratio is
            #an algebraic integer, which equal norms do NOT imply: for alpha = -2
            #both 1+x and 1-x have norm 3 but their ratio is (-1+2x)/3.  When the
            #ratio is not integral this falls through to the general Euclidean
            #step below, which is what the norm-Euclidean property provides.
            q = matric[0][0]*(matric[0][1])^(-1)
            matric = matric*Afunc(-q)
            relatorlist.insert(0,[f(q)[0],f(q)[1],0])
            continue#So far the 'very close to stopping ones'
        elif  norms(matric[0][0]) < norms(matric[0][1]):
            relatorlist.insert(0,[0,0,3])
            matric = matric*S
        else:
            #print(matric)
            z = Euclidian_Quotient(matric[0][0],matric[0][1])
            #print(z)
            vecz = f(z)
            matric = matric*Afunc(-z)
            #The loop right-multiplies until the identity, so M = gn^-1 ... g1^-1
            #and each recorded syllable must be the INVERSE of the multiplier.
            #The multiplier here is A^-z, so record A^z.
            relatorlist.insert(0,[vecz[0],vecz[1],0])
    relatorlist = concatenate_loop(relatorlist)
    assert(word_to_vector(copy.deepcopy(relatorlist)) == Target_matrix)
    return [GLFactor,relatorlist]

def Commuting_square_relations():
    #At a prime-power level with exponent >= 3, an element
    #A^x S A^(1+omega) S^-1 A^y with x,y NOT coprime to the level can be written
    #two ways: stepping through the real part first, or through the imaginary
    #part first.  Same element, two different words.  PriPow_relation always
    #takes the real-first path, so the difference of the two words is a relator
    #that Reidemeister-Schreier never emits -- which is why the relation lattice
    #came out with too small a rank and the Hecke operators failed to commute.
    #At exponent <= 2 every zero-divisor relator has the same front and tail
    #(A^x S A S^-1 A^-x), so the ordering cannot matter and these relations come
    #out identically zero -- verified, so adding them changes nothing there.
    #Rank restored: e=3 22 -> 24, e=4 64 -> 72, for d=-2 and d=-11 alike.
    pieces = Seperate_PriPow_generators(Real_Gens)
    if pieces == []:
        return []
    real_gens, im_gens = pieces[0], pieces[1]
    def start_of(g):
        return g[0][0]
    def end_of(g):
        return [-g[0][2][0], -g[0][2][1], 1]
    squares = []
    for r1 in real_gens:
        s = start_of(r1)
        e1 = end_of(r1)
        for m1 in im_gens:
            if start_of(m1) != e1:
                continue
            target = end_of(m1)
            for m2 in im_gens:
                if start_of(m2) != s:
                    continue
                e2 = end_of(m2)
                for r2 in real_gens:
                    if start_of(r2) != e2:
                        continue
                    if end_of(r2) == target:
                        squares.append([[r1[1],1],[m1[1],1],[m2[1],-1],[r2[1],-1]])
    return squares

def All_Relations(PSL, GL):

    num = len(Real_Gens)
    if alpha == -1:
        half_unit_ord = 2
    elif alpha == -3:
        half_unit_ord = 3
    else:
        half_unit_ord = 1

    The_Unit_rel = Unit_Relation()
    All_Rels = [The_Unit_rel]
    All_Rels += Fin_Oder_Inv_gen()
    for rel in Real_Rels:
        All_Rels.append(Rewrite_Relation(rel))
    #the real-first / imaginary-first ambiguity at exponent >= 3 (zero at e <= 2)
    All_Rels += Commuting_square_relations()
    if GL == True:
        All_Rels.append([[num,2*half_unit_ord]])
        GLRels = GL_Relator_Strings()
        for rel in GLRels:
            All_Rels.append(Rewrite_Relation_GL(rel))
        if PSL == True:
            All_Rels.append([[0, half_unit_ord]])
            All_Rels.append([[0,1],[num,2]])
        else:
            All_Rels.append([[0,2*half_unit_ord]])
    elif PSL == True:
        All_Rels.append([[num,1]])
        All_Rels.append([[0,half_unit_ord]])
    else:
        All_Rels.append([[num,1]])
        All_Rels.append([[0,2*half_unit_ord]])
    return All_Rels

def relatormatrix(rewritten):
    relations = rewritten
    gl_rels = []
    matrixfromvec = []
    n = len(Real_Gens) + 1
    for i in relations:
        vec = [0] * n
        #print(i, "WHOOP")
        for j in range(len(i)):
            vec[i[j][0]] += i[j][1]
        if vec not in matrixfromvec:
            matrixfromvec.append(vec)
    Returnmatrix = Matrix(matrixfromvec)
    return Returnmatrix    
def Nottransfer_takeout(matrx, columnindex):
    i = columnindex
    num = GCD(matrx.column(i))
    m = modulus(matrx.column(i))
    j = smallest(matrx,i)
    while m > num^2:
        j = smallest(matrx,i)
        n1 = matrx.nrows()
        for k in range(n1):
            if k != j:
                q1 = floor(matrx[k][i]/matrx[j][i])
                q2 = q1+1
                a1 = matrx[k][i] - q1*matrx[j][i]
                a2 = matrx[k][i] - q2*matrx[j][i]
                if abs(a2) < abs(a1):
                    matrx[k] = matrx[k] - q2*matrx[j]
                else:
                    matrx[k] = matrx[k] - q1*matrx[j]
        m = modulus(matrx.column(i))
    matrx = matrx.delete_rows([j])
    matrx = matrx.delete_columns([i])
    #matrx = row_column_remover(matrx)
    return matrx

def Nottransfer_final(matrx):
    n1 = matrx.ncols()
    for i in reversed(range(n1)):
        if GCD(matrx.column(i)) == 1:
            matrx = Nottransfer_takeout(matrx, i)
            #print(matrx)
    for i in matrx.columns():
        assert(GCD(i) != 1)
    return matrx

def nearest_quotient(a, b):
    #Rounds a/b to the nearest integer, breaking ties downwards.
    #Same "try floor, then floor+1" trick already used in Nottransfer_takeout.
    q1 = floor(a/b)
    q2 = q1 + 1
    if abs(a - q2*b) < abs(a - q1*b):
        return q2
    return q1

def smallest_among(matrx, colj, rowlist):
    #Like smallest(), but only looks at the rows in rowlist.
    #Returns the row with the smallest non-zero absolute value in column colj,
    #ties broken by the row of smallest modulus.
    index = -1
    temp = 0
    m = 0
    for i in rowlist:
        entry = matrx[i][colj]
        if entry == 0:
            continue
        if index == -1 or abs(entry) < temp:
            temp = abs(entry)
            index = i
            m = modulus(list(matrx[i]))
        elif abs(entry) == temp and modulus(list(matrx[i])) < m:
            index = i
            m = modulus(list(matrx[i]))
    return index

def Column_Gcd_Reduce(matrx, verbose=False, max_sweeps=50):
    #Sweeps the columns of an integral matrix from right to left.  For each
    #column the gcd of the entries is computed and, by integer row operations,
    #the column is reduced until a single row carries that gcd, made positive.
    #
    #Reducing a column re-introduces entries into the columns to its right, so
    #one sweep is not enough and the sweep is repeated.  Row operations alone
    #cannot always reach one entry per column (that is Smith normal form, which
    #needs column operations too), so the repetition stops at the best state it
    #can find rather than at a guaranteed diagonal.
    #
    #A row is frozen only when the gcd of its column came out as 1, i.e. only
    #when that row and column are going to be taken out at the end.  Freezing
    #every pivot instead would shrink the pool of rows available to the columns
    #further left and make gcds that ought to reach 1 come out too big.
    #
    #The matrix is kept at its original size throughout and the frozen rows and
    #columns are deleted only at the very end, together with every all-zero row.
    #
    #Returns (reduced_matrix, taken_rows).  taken_rows still carries all of the
    #original columns, so the position of its 1 records which column it killed.
    mat = matrix(ZZ, matrx)
    nr = mat.nrows()
    nc = mat.ncols()
    take_rows = []
    take_cols = []

    best = copy.deepcopy(mat)
    best_score = sum([modulus(list(mat[i])) for i in range(nr)])
    for sweep in range(max_sweeps):
        before = copy.deepcopy(mat)
        taken_before = len(take_rows)
        for j in reversed(range(nc)):
            if j in take_cols:
                continue
            free = [r for r in range(nr) if r not in take_rows]
            live = [r for r in free if mat[r][j] != 0]
            if len(live) == 0:
                continue
            g = GCD([mat[r][j] for r in live])
            while len(live) > 1:
                p = smallest_among(mat, j, live)
                for r in live:
                    if r != p:
                        mat[r] = mat[r] - nearest_quotient(mat[r][j], mat[p][j])*mat[p]
                live = [r for r in free if mat[r][j] != 0]
            p = live[0]
            assert abs(mat[p][j]) == g
            #Force the surviving entry to be +gcd.  Without this the sign of a
            #pivot can flip from one sweep to the next and the loop below never
            #sees two equal sweeps in a row.
            if mat[p][j] < 0:
                mat[p] = -mat[p]
            if verbose:
                print('sweep', sweep, 'column', j, 'gcd', g, 'pivot row', p)
            if g != 1:
                continue
            #mat[p] is zero in every column taken out so far, so adding multiples
            #of it to those frozen rows clears column j without disturbing them.
            for r in take_rows:
                if mat[r][j] != 0:
                    mat[r] = mat[r] - nearest_quotient(mat[r][j], mat[p][j])*mat[p]
            take_rows.append(p)
            take_cols.append(j)
        #A sweep can undo part of an earlier one, so the loop is not guaranteed
        #to settle by itself.  Carry on only while the entries genuinely get
        #smaller (or a row was just taken out), and hand back the best state
        #seen.  That makes the total strictly decrease, so this always stops.
        if mat == before:
            break
        if len(take_rows) > taken_before:
            best = copy.deepcopy(mat)
            best_score = sum([modulus(list(mat[i])) for i in range(nr)])
            continue
        score = sum([modulus(list(mat[i])) for i in range(nr)])
        if score < best_score:
            best = copy.deepcopy(mat)
            best_score = score
        else:
            mat = best
            break

    for j, r in zip(take_cols, take_rows):
        col = mat.column(j)
        assert GCD(col) == 1
        assert [i for i in range(nr) if col[i] != 0] == [r]
        if mat[r][j] == -1:
            mat[r] = -mat[r]

    taken = mat.matrix_from_rows(take_rows)
    #Keep the eliminated column indices: they say which generators were made
    #redundant, which Reduce_vector needs to push a vector onto the survivors.
    Taken_columns = list(take_cols)
    reduced = mat.delete_rows(take_rows).delete_columns(take_cols)
    #Drop the rows that carry nothing.  An all-zero row of the input is never
    #touched above (it is never "live" in any column), so it is still zero here
    #and gets dropped; rows that turned zero during the reduction go too.
    reduced = reduced.delete_rows([i for i in range(reduced.nrows())
                                   if reduced[i].is_zero()])
    return reduced, taken, Taken_columns

def Return_Deg1_Prime(prime_number):
    factorlist = list(factor(O(prime_number)))
    #print(factorlist)
    if len(factorlist) == 2 or factorlist[0][1] == 2:
        return True
    else:
        return False



def Reduce_vector(vec, taken, taken_columns):
    #Rewrite a relation vector so that it is supported away from the generators
    #Column_Gcd_Reduce eliminated.  Each taken row carries a +1 in its own taken
    #column and 0 in every other taken column, so subtracting vec[j]*row zeroes
    #column j without disturbing the ones already done: one pass, order-free.
    #The result represents the same class in the abelianisation.
    #The coefficient ring is read off the taken rows, so this serves the integral
    #reduction (ZZ) and the field reductions (QQ or GF(p)) alike.
    R = taken.base_ring()
    v = vector(R, vec)
    for k in range(len(taken_columns)):
        j = taken_columns[k]
        if v[j] != 0:
            v = v - v[j]*vector(R, taken[k])
    return v

def Field_Gcd_Reduce(matrx, K):
    #The field analogue of Column_Gcd_Reduce.  Over a field every nonzero column
    #has a pivot, so the gcd sweep collapses to reduced row echelon form -- which
    #is Sage's echelon_form, in C, instead of a Python sweep.  Reducing into K
    #before sweeping also keeps the entries small, which is where the time goes
    #for large index.  Sweeping columns in reverse is the same as echelonising
    #the matrix with its columns reversed.
    mat = matrix(K, matrx)
    nc = mat.ncols()
    backwards = list(reversed(range(nc)))
    echelon = mat.matrix_from_columns(backwards).echelon_form()
    pivots = list(echelon.pivots())
    taken_columns = [nc-1-j for j in pivots]
    taken = echelon.matrix_from_rows(range(len(pivots))).matrix_from_columns(backwards)
    #every non-pivot row of an echelon form vanishes, so no relations survive:
    #the abelianisation tensored with K is free on the surviving generators
    reduced = matrix(K, 0, nc - len(pivots))
    return reduced, taken, taken_columns

def Surviving_generators(ncols, taken_columns):
    #The generator indices that survive the reduction, in order.
    return [j for j in range(ncols) if j not in taken_columns]

def Conjugate_elt(z):
    #Complex conjugation on O: x -> 1-x when alpha = 1 (mod 4), else x -> -x.
    v = f(O(z))
    if alpha % 4 == 1:
        return O(v[0]) + O(v[1]) - O(v[1])*O(x)
    return O(v[0]) - O(v[1])*O(x)

def Rational_below(pi):
    #The rational prime lying under the prime pi.
    return ZZ(O(pi).absolute_norm()).factor()[0][0]

def Prime_residues(pi):
    #A complete residue system for O/pi with pi prime: the rational integers
    #0..p-1 when pi has degree one, and a+b*x with 0 <= a,b < p when pi is
    #inert.  Same shape as the transversal description on p.14 of the thesis.
    p = Rational_below(pi)
    if O(pi).absolute_norm() == p:
        return [O(i) for i in range(p)]
    return [O(i) + O(j)*O(x) for i in range(p) for j in range(p)]

def Hecke_cosets(pi):
    #Right-coset representatives for Gamma_pi \ Gamma_0(level), as in Lemma 3.2.
    #Here Gamma_pi = Gamma_0(level) intersect delta_pi Gamma_0(level) delta_pi^-1,
    #and conjugating shows that is exactly Gamma_0(pi*level), so the index is
    #N(pi) when pi divides the level and N(pi)+1 otherwise.
    #
    #The A^{level*s} family is common to all three cases and is the whole
    #transversal when pi divides the level.  Otherwise one more representative
    #A^e S A^y is needed, with
    #   e = the rational prime below pi                     (general case),
    #   e = the generator of the conjugate ideal            (pi conjugate to
    #                                                        the level's prime),
    #and y the inverse of e modulo the level, so that 1 - e*y = 0 mod level.
    I   = O.fractional_ideal(level)
    Ipi = O.fractional_ideal(O(pi))
    reps = [Afunc(level*s) for s in Prime_residues(pi)]
    if I + Ipi == Ipi:
        return reps
    conj_base = Conjugate_elt(prime_base)
    if Ipi == O.fractional_ideal(O(conj_base)):
        e = O(conj_base)
    else:
        e = O(Rational_below(pi))
    y = O(e).inverse_mod(I)
    reps.append(Afunc(e)*S*Afunc(y))
    return reps

Residue_field_cache = {}
Divides_level_cache = {}
Transversal_cache = {}

def Hecke_transversal(pi):
    #the coset representatives together with their index map, built once per
    #(level, pi).  Hecke_index_map costs one coset key per representative, so
    #rebuilding it per generator (or per harvested column) was pure waste.
    key = (str(level), str(O(pi)))
    if key not in Transversal_cache:
        reps = Hecke_cosets(pi)
        Transversal_cache[key] = (reps, Hecke_index_map(pi))
    return Transversal_cache[key]

def Residue_field(pi):
    #O/pi as a genuine finite field, with a fast coercion.  The Hecke operators
    #are indexed by PRIME elements, so this always exists, and arithmetic in it
    #replaces the ideal operations (reduce / inverse_mod / membership) that
    #dominated the profile of Hecke_coset_index.
    key = str(O(pi))
    if key not in Residue_field_cache:
        Ipi = O.fractional_ideal(O(pi))
        assert Ipi.is_prime(), '%s does not generate a PRIME ideal -- Hecke operators are indexed by primes, and the coset key is computed in the residue field O/pi.  Factor it first (e.g. via Good_primes).' % O(pi)
        Residue_field_cache[key] = Ipi.residue_field()
    return Residue_field_cache[key]

def Divides_level(pi):
    key = (str(level), str(O(pi)))
    if key not in Divides_level_cache:
        Ipi = O.fractional_ideal(O(pi))
        Divides_level_cache[key] = (O.fractional_ideal(level) + Ipi == Ipi)
    return Divides_level_cache[key]

def Hecke_coset_index(M, pi, index_map=None, kf=None, divides=None):
    #Canonical key for the coset of M, computed in the residue field O/pi: a
    #point of P^1(O/pi) is determined by the ratio of the bottom row, and
    #d invertible mod pi is just kf(d) != 0.  Same keys as the ideal-arithmetic
    #version, but without constructing ideals per call.
    #kf and divides are hoisted out of the caller's loop: looking them up here
    #rebuilds the cache key with str(O(pi)) on EVERY coset, which profiling showed
    #costs more than the arithmetic it guards.
    if kf is None:
        kf = Residue_field(pi)
    if divides is None:
        divides = Divides_level(pi)
    #coerce the bottom row STRAIGHT into O/pi.  Wrapping in O(...) first was the
    #top line of the profile (order._element_constructor_); the residue field
    #accepts number-field elements that are integral at pi, which these are.
    c = M[1][0]; d = M[1][1]
    dbar = kf(d)
    if divides:
        key = ('L', kf(c/level)/dbar)
    elif dbar != 0:
        key = ('A', kf(c)/dbar)
    else:
        key = ('S', dbar/kf(c))
    if index_map is None:
        return key
    return index_map[key]

def Hecke_index_map(pi):
    #key -> position in Hecke_cosets(pi); build once and reuse.
    kf = Residue_field(pi)
    divides = Divides_level(pi)
    return dict((Hecke_coset_index(M, pi, None, kf, divides), j) for j, M in enumerate(Hecke_cosets(pi)))

def Unit_generator():
    #the unit u used by the unit Hecke operator: x for alpha = -1,-3 (i and a
    #primitive 6th root of unity), -1 otherwise.
    if alpha in (-1,-3):
        return x
    return -1

def Relation_to_vector(rel, n):
    #a rewritten word [[gen,exponent],...] as its exponent vector
    v = [0]*n
    for term in rel:
        v[term[0]] += term[1]
    return vector(ZZ, v)

def Element_to_vector(M, n):
    #matrix in Gamma_0(level) -> word (Euclidian algorithm) -> rewritten word
    #in the subgroup generators -> exponent vector
    word = Euclidian_algorithm(copy.deepcopy(M))[1]
    return Relation_to_vector(Rewrite_Relation_GL(copy.deepcopy(word)), n)

def Reduction_ring(coefficients):
    #-1 = no reduction (the integral abelianisation), 0 = the free part over QQ,
    #a rational prime p = the mod p part over GF(p)
    if coefficients == -1:
        return ZZ
    if coefficients == 0:
        return QQ
    return GF(coefficients)

def Abelianisation_reduction(PSL=False, GL=False, coefficients=-1):
    #The relation lattice together with its reduction, so that every Hecke
    #operator is expressed in one and the same basis of surviving generators.
    #Compute this once per level and reuse it -- the sweep is the expensive part.
    #`coefficients` picks the coefficient ring: -1 keeps everything integral and
    #gives the full abelianisation; 0 works over QQ (the free part); a rational
    #prime p works over GF(p).  Reducing at the START of the sweep is what makes
    #the large-index cases cheap.  The bundle is stamped with the level and the
    #flags it was built from, so a stale one cannot be used by accident.
    assert coefficients in [-1, 0] or is_prime(coefficients), 'coefficients must be -1, 0 or a rational prime'
    Rel = matrix(ZZ, relatormatrix(All_Relations(PSL, GL)))
    n = Rel.ncols()
    if coefficients == -1:
        reduced, taken, taken_columns = Column_Gcd_Reduce(Rel)
    else:
        reduced, taken, taken_columns = Field_Gcd_Reduce(Rel, Reduction_ring(coefficients))
    surviving = Surviving_generators(n, taken_columns)
    return [reduced, taken, taken_columns, surviving, n, PSL, GL, level, coefficients]

def Current_reduction(GL=False, reduction=None, coefficients=None):
    #resolve which reduction to use: an explicitly supplied one, else the global
    #Fixed_Reduction, else a fresh one.  coefficients=None means "whatever the
    #resolved reduction already uses", so setting a mod p global makes everything
    #downstream work mod p without further argument passing.
    if reduction is None:
        reduction = globals().get('Fixed_Reduction')
    if reduction is None:
        return Abelianisation_reduction(False, GL, -1 if coefficients is None else coefficients)
    #a level mismatch means the global is stale -- that is an error.  A flag or
    #coefficient mismatch just means a different operator was asked for.
    assert reduction[7] == level, 'Fixed_Reduction was built for level %s, not %s -- recompute it' % (reduction[7], level)
    if coefficients is None:
        coefficients = reduction[8]
    if reduction[6] != GL or reduction[8] != coefficients:
        return Abelianisation_reduction(reduction[5], GL, coefficients)
    return reduction

def All_Mat_Gens(GL=False):
    #the matrix generators; in the GL case the extra generator diag(1,u) sits at
    #the last index, matching the extra column of the GL relation matrix
    if GL:
        return Real_Mat_Gens + [Matrix([[1,0],[0,Unit_generator()]])]
    return Real_Mat_Gens

def GElement_to_vector(M, n):
    #as Element_to_vector, but for GL: split off diag(1,u)^k where det(M) = u^k,
    #decompose the determinant-1 part, then record k on the last generator
    u = Unit_generator()
    determinant = M.determinant()
    k = None
    for t in range(2*half_unit_ord):
        if determinant == u^t:
            k = t
            break
    assert k is not None, 'determinant is not a power of the unit'
    v = list(Element_to_vector(M*(Matrix([[1,0],[0,u]])^(-k)), n))
    v[n-1] += k
    return vector(ZZ, v)

def Unit_Hecke_matrix(GL=False, reduction=None, coefficients=None):
    #The unit Hecke operator T_u.  Because u is a unit, conjugation by
    #delta_u = diag(1,u) preserves Gamma_0(level), so Gamma_u = Gamma: the
    #transversal is a single coset and the transfer is trivial.  T_u therefore
    #collapses to plain conjugation g -> delta_u^-1 g delta_u, in the same
    #direction as T_pi, followed by the usual decomposition and reduction.
    #Sanity check: with GL=True the conjugation is by a group element, hence
    #inner, so T_u must come out as the identity.  With GL=False it is a
    #genuine involution (T_u^2 = id, Lemma 3.4) but not the identity.
    reduction = Current_reduction(GL, reduction, coefficients)
    taken, taken_columns, surviving, n = reduction[1], reduction[2], reduction[3], reduction[4]
    generators = All_Mat_Gens(GL)
    to_vector = GElement_to_vector if GL else Element_to_vector
    delta = Matrix([[1,0],[0,Unit_generator()]])
    delta_inv = delta.inverse()
    columns = []
    for k in surviving:
        conjugated = delta_inv*generators[k]*delta
        v = Reduce_vector(to_vector(conjugated, n), taken, taken_columns)
        columns.append([v[j] for j in surviving])
    return Matrix(Reduction_ring(reduction[8]), columns).transpose()

def Hecke_operator_vector(g, pi, GL=False, reduction=None, reps=None, index_map=None):
    #One column of T_pi, before reduction: the transfer of g from Gamma to
    #Gamma_pi = Gamma_0(pi*level), conjugated into Gamma by delta_pi = diag(1,pi)
    #and summed over the transversal.  This is the composite t, then delta, then
    #i of diagram (8): for each coset rep h_i, the element h_i g h_j^-1 lies in
    #Gamma_pi, and delta_pi^-1 (h_i g h_j^-1) delta_pi lands back in Gamma.
    #reps/index_map are passed in by Hecke_matrix so the transversal is built once
    #per prime instead of once per generator.
    n = reduction[4]
    to_vector = GElement_to_vector if GL else Element_to_vector
    if reps is None or index_map is None:
        reps, index_map = Hecke_transversal(pi)
    delta = Matrix([[1,0],[0,O(pi)]])
    delta_inv = delta.inverse()
    total = vector(ZZ, [0]*n)
    kf = Residue_field(pi)
    divides = Divides_level(pi)
    for i in range(len(reps)):
        product = reps[i]*g
        j = Hecke_coset_index(product, pi, index_map, kf, divides)
        element = product*(reps[j].inverse())
        conjugated = delta_inv*element*delta
        if i == 0:
            #delta^-1 M delta is integral for every M in Gamma_pi (pi divides the
            #lower-left entry), so this is a theorem, not a runtime condition.
            #Checking it once per prime catches a systematic error without paying
            #for a membership test on every coset -- that test alone was a third
            #of the running time.
            assert all(entry in O for entry in conjugated.list()), 'conjugated element left the ring'
        total = total + to_vector(conjugated, n)
    return total

Hecke_matrix_cache = {}
Hecke_vector_cache = {}

def Clear_Hecke_cache():
    #call this after changing level or regenerating the group data
    global Hecke_matrix_cache, Hecke_vector_cache, Residue_field_cache, Divides_level_cache, Transversal_cache
    Hecke_matrix_cache = {}
    Hecke_vector_cache = {}
    Residue_field_cache = {}
    Divides_level_cache = {}
    Transversal_cache = {}

def Hecke_matrix(pi, GL=False, reduction=None, coefficients=None):
    #The Hecke operator T_pi on the abelianisation, in the basis of surviving
    #generators.  One column per surviving generator.  Cached, because the
    #decomposition routines ask for the same operator repeatedly and each rebuild
    #costs (N(pi)+1) * (#generators) Euclidean decompositions.
    reduction = Current_reduction(GL, reduction, coefficients)
    key = (str(level), reduction[5], GL, reduction[8], str(O(pi)))
    if key in Hecke_matrix_cache:
        return Hecke_matrix_cache[key]
    taken, taken_columns, surviving = reduction[1], reduction[2], reduction[3]
    generators = All_Mat_Gens(GL)
    reps, index_map = Hecke_transversal(pi)
    columns = []
    for k in surviving:
        #the raw transfer vector is INTEGRAL and does not depend on the
        #coefficient ring -- only the reduction that follows does.  Caching it
        #makes running the same level over QQ and over each torsion prime cost
        #one set of Euclidean decompositions instead of one set per ring.
        vkey = (str(level), GL, str(O(pi)), k)
        if vkey in Hecke_vector_cache:
            total = Hecke_vector_cache[vkey]
        else:
            total = Hecke_operator_vector(generators[k], pi, GL, reduction, reps, index_map)
            Hecke_vector_cache[vkey] = total
        v = Reduce_vector(total, taken, taken_columns)
        columns.append([v[j] for j in surviving])
    T = Matrix(Reduction_ring(reduction[8]), columns).transpose()
    Hecke_matrix_cache[key] = T
    return T

def Elementary_divisors(PSL=False, GL=False, reduction=None):
    #The elementary divisors of the abelianisation.  Taken from the INTEGRAL
    #reduction: Column_Gcd_Reduce only ever removes a column with a unit pivot
    #together with its row, which is a unimodular change of basis, so the
    #cokernel is unchanged and the Smith form of the small reduced matrix is the
    #Smith form of the whole thing -- much cheaper than smithing the full matrix.
    reduction = Current_reduction(GL, reduction, -1)
    reduced = reduction[0]
    if reduced.nrows() == 0 or reduced.ncols() == 0:
        return []
    D = reduced.smith_form(False)
    return [D[i][i] for i in range(min(D.nrows(), D.ncols())) if D[i][i] not in [0,1]]

def Torsion_primes(PSL=False, GL=False, reduction=None):
    #Exactly the primes at which the mod p reduction sees more than the free
    #part does, i.e. the primes worth passing as `coefficients`.  For any other
    #prime the mod p dimension equals the free rank and nothing new appears.
    divisors = Elementary_divisors(PSL, GL, reduction)
    return sorted(set([q for d in divisors for (q, e) in factor(d)]))

def Relevant_coefficients(PSL=False, GL=False, reduction=None):
    #The full list of coefficient settings worth running: 0 for the free part,
    #then one entry per torsion prime.
    return [0] + Torsion_primes(PSL, GL, reduction)

def Phi_primes():
    #The primes dividing the order of (O/level)^x, the target of the map phi.
    #This is ONLY a heuristic for where an Eisenstein class may appear.  It is
    #NOT a classification of the prime:
    #  - at such a prime there are still non-trivial NON-Eisenstein eigenclasses;
    #  - the Eisenstein contribution need not sit in the torsion at all, it can
    #    come from a free part that is not the trace class.
    #So never label a prime "Eisenstein"; decompose the space and classify the
    #eigenCLASSES with Eisenstein_subspace / Trace_subspace below.
    N = ZZ(O(level).absolute_norm())
    p_norm = ZZ(O(prime_base).absolute_norm())
    order = (N//p_norm)*(p_norm - 1)
    return [q for (q, e) in factor(order)]

def Hecke_eigenspace(eigenvalue, pis, GL=False, reduction=None, coefficients=None):
    #The simultaneous eigenspace on which T_pi acts by the scalar eigenvalue(pi),
    #for every pi in the list: the intersection of the right kernels of
    #T_pi - eigenvalue(pi).  Field coefficients only (0 or a prime), since this
    #is linear algebra over a field.  Hecke_matrix acts on column vectors, so the
    #right kernel is the correct one.
    reduction = Current_reduction(GL, reduction, coefficients)
    assert reduction[8] != -1, 'Hecke_eigenspace needs field coefficients: 0 for the free part, or a prime'
    K = Reduction_ring(reduction[8])
    r = len(reduction[3])
    W = VectorSpace(K, r)
    for pi in pis:
        T = Hecke_matrix(pi, GL, reduction)
        W = W.intersection((T - K(eigenvalue(pi))*identity_matrix(K, r)).right_kernel())
        if W.dimension() == 0:
            break
    return W

def Eisenstein_subspace(pis, GL=False, reduction=None, coefficients=None):
    #the classes with a_pi = N(pi)+1 for every pi tested, wherever they live
    return Hecke_eigenspace(lambda pi: ZZ(O(pi).absolute_norm())+1, pis, GL, reduction, coefficients)

def Trace_subspace(pis, GL=False, reduction=None, coefficients=None):
    #the classes with a_pi = Tr(pi) for every pi tested
    return Hecke_eigenspace(lambda pi: ZZ(O(pi).trace()), pis, GL, reduction, coefficients)

def Good_primes(bound=20, include_level=True):
    #Prime elements to test with.  Primes DIVIDING the level are included by
    #default: U_pi at the level is what separates a class from its Galois
    #conjugate -- its characteristic polynomial there is the minimal polynomial
    #of pi -- so omitting them leaves the trace class merged as a 2-dimensional
    #piece instead of splitting into two conjugate lines.
    out = []
    for p in primes(bound):
        for (pi_, e_) in factor(O(p)):
            pi = O(pi_)
            divides = (O.fractional_ideal(level) + O.fractional_ideal(pi) == O.fractional_ideal(pi))
            if divides and not include_level:
                continue
            out.append(pi)
    return out

def Restrict_operator(T, W):
    #The matrix of T acting on the T-stable subspace W, in W's own basis.
    #Hecke_matrix acts on COLUMN vectors, so the basis vectors go in as columns.
    B = W.basis_matrix().transpose()
    return B.solve_right(T*B)

def Lift_subspace(W, U):
    #U is given in W's coordinates; return it as a subspace of the ambient space
    B = W.basis_matrix().transpose()
    if U.dimension() == 0:
        return (B*matrix(B.base_ring(), B.ncols(), 0)).column_space()
    return (B*U.basis_matrix().transpose()).column_space()

def Refine_decomposition(pieces, pi, GL=False, reduction=None):
    #Split each current piece along the primary components of ONE more operator.
    #Primary (ker g(T)^m) rather than ker g(T), so a non-semisimple action loses
    #nothing.
    T = Hecke_matrix(pi, GL, reduction)
    refined = []
    for W in pieces:
        if W.dimension() <= 1:
            refined.append(W)
            continue
        M = Restrict_operator(T, W)
        for (g, m) in M.charpoly().factor():
            U = (g(M)^m).right_kernel()
            if U.dimension() > 0:
                refined.append(Lift_subspace(W, U))
    return refined

def Hecke_decomposition(pis, GL=False, reduction=None, coefficients=None):
    #The finest decomposition the given primes can see; each piece carries a
    #single eigenvalue system (possibly a Galois orbit of one).
    reduction = Current_reduction(GL, reduction, coefficients)
    K = Reduction_ring(reduction[8])
    assert reduction[8] != -1, 'needs field coefficients: 0 for the free part, or a prime'
    pieces = [VectorSpace(K, len(reduction[3]))]
    for pi in pis:
        pieces = Refine_decomposition(pieces, pi, GL, reduction)
    return pieces

def Generic_combination(pis, W, GL=False, reduction=None, tries=12):
    #A linear combination of the Hecke operators whose eigenspaces on W are the
    #simultaneous ones.  Generic coefficients work; we just try a few small
    #deterministic choices and keep the first that is as split as possible, i.e.
    #whose characteristic polynomial is squarefree on W.
    reduction = Current_reduction(GL, reduction, None)
    mats = [Restrict_operator(Hecke_matrix(pi, GL, reduction), W) for pi in pis]
    K = mats[0].base_ring()
    best = mats[0]
    for t in range(tries):
        coeffs = [K(1 + ((t+1)*(i+1)) % 7) for i in range(len(mats))]
        R = sum([coeffs[i]*mats[i] for i in range(len(mats))])
        if R.charpoly().is_squarefree():
            return R, mats
        best = R
    return best, mats

def Eigenvalue_system(pis, W, GL=False, reduction=None):
    #The eigenvalue system on a primary piece W, over an extension of the base
    #field, returned as [L, system, eigenvector (ambient), eigenspace dim].
    #
    #Two things have to be handled at once, and getting either wrong shows up as
    #a spurious "not a single eigenvalue system":
    #  * NON-SEMISIMPLICITY: T_pi = a_pi + nilpotent, so a vector killed by the
    #    nilpotent part of a generic combination R need not be killed by that of
    #    T_pi.  Reading a_pi off such a vector fails.
    #  * GALOIS ORBITS: a piece irreducible over K splits over the residue field
    #    L into [L:K] conjugate eigenlines, so T_pi has SEVERAL roots on W over L.
    #    Demanding a single root would reject a perfectly good orbit.
    #Both are solved by passing to the generalised eigenspace of R for ONE chosen
    #root: that fixes an embedding (killing the orbit ambiguity) and contains the
    #whole Jordan block (killing the nilpotent ambiguity).  On it every T_pi has
    #exactly one eigenvalue.
    R, mats = Generic_combination(pis, W, GL, reduction)
    K = R.base_ring()
    g = R.charpoly().factor()[0][0]
    if g.degree() == 1:
        L = K
        root = -g[0]/g[1]
    elif K == QQ:
        L = NumberField(g, 'a')
        root = L.gen()
    else:
        L = K.extension(g.degree(), 'a')
        root = g.change_ring(L).roots()[0][0]
    r = R.nrows()
    RL = R.change_ring(L)
    G = ((RL - root*identity_matrix(L, r))^r).right_kernel()
    assert G.dimension() >= 1, 'empty generalised eigenspace'
    BG = G.basis_matrix().transpose()
    system = []
    restricted = []
    for M in mats:
        ML = M.change_ring(L)
        MG = BG.solve_right(ML*BG)          #restrict to the generalised eigenspace
        roots = MG.charpoly().roots(L)
        assert len(roots) == 1, 'this piece carries more than one eigenvalue system -- add more primes'
        system.append(roots[0][0])
        restricted.append(MG)
    #common eigenvector inside the generalised eigenspace
    V = VectorSpace(L, G.dimension())
    for i in range(len(restricted)):
        V = V.intersection((restricted[i] - system[i]*identity_matrix(L, G.dimension())).right_kernel())
    assert V.dimension() >= 1, 'no common eigenvector -- the operators do not commute here'
    v = BG*V.basis()[0]
    eigenvector = W.basis_matrix().transpose().change_ring(L)*v
    if L != K and K == QQ:
        try:
            Lopt, from_opt, to_opt = L.optimized_representation()
            system = [to_opt(a) for a in system]
            eigenvector = vector([to_opt(c) for c in eigenvector])
            L = Lopt
        except Exception:
            pass
    return [L, system, eigenvector, V.dimension()]

def Hecke_eigensystems(pis, GL=False, reduction=None, coefficients=None):
    #Full report, one entry per primary piece:
    #  [generalised dim, eigenspace dim, residue field, {pi: a_pi}, eigenvector]
    #generalised dim is the dimension of the primary component the eigenvector
    #lives in.  Note that gen.dim > eig.dim does NOT by itself mean the action is
    #non-semisimple: a piece irreducible over K splits over the residue field L
    #into [L:K] conjugate eigenlines, so the semisimple case is
    #    gen.dim == [L:K] * eig.dim
    #and only a strict excess over that is a genuine Jordan block.
    reduction = Current_reduction(GL, reduction, coefficients)
    out = []
    for W in Hecke_decomposition(pis, GL, reduction):
        L, system, eigenvector, eigdim = Eigenvalue_system(pis, W, GL, reduction)
        out.append([W.dimension(), eigdim, L, dict(zip([str(p) for p in pis], system)), eigenvector])
    return out

def Primes_by_norm(bound=40, include_level=True, max_norm=None):
    #Primes sorted by NORM, cheapest first.  Building T_pi costs (N(pi)+1)
    #element decompositions, so norm IS the cost.  Note an inert rational prime p
    #has N = p^2, so it is p times dearer than a degree-one prime above p: for
    #bulk eigenvalue data prefer the degree-one primes.
    out = Good_primes(bound, include_level)
    if max_norm is not None:
        out = [pi for pi in out if ZZ(O(pi).absolute_norm()) <= max_norm]
    return sorted(out, key=lambda pi: ZZ(O(pi).absolute_norm()))

def Verify_right_eigenvector(entry, pis, GL=False, reduction=None):
    #Check directly that the returned vector is a RIGHT eigenvector: T_pi v = a_pi v.
    reduction = Current_reduction(GL, reduction, None)
    L, v = entry[2], entry[4]
    for pi in pis:
        T = Hecke_matrix(pi, GL, reduction).change_ring(L)
        a = L(entry[3][str(pi)])
        if T*v != a*v:
            return False
    return True

def Certified_line(entry):
    #entry is one row of Hecke_eigensystems.  The eigenspace being 1-dimensional
    #over the residue field is a CERTIFICATE: for W = intersection of the
    #ker(T_pi - a_pi) over the tested primes, any Hecke operator T_rho commutes
    #with those T_pi, hence preserves W, hence acts on a 1-dimensional W by a
    #scalar.  So no further eigenvector check is needed for any new prime -- the
    #eigenvalues can simply be harvested.  This holds even when the surrounding
    #generalised eigenspace is a Jordan block, since T_rho preserves the honest
    #kernel too.
    return entry[1] == 1

def Hecke_column(pi, j, GL=False, reduction=None):
    #Column j of T_pi (j indexing the surviving generators), without building the
    #rest of the matrix.  One column costs (N(pi)+1) element decompositions
    #instead of r*(N(pi)+1) for the full operator.
    reduction = Current_reduction(GL, reduction, None)
    taken, taken_columns, surviving = reduction[1], reduction[2], reduction[3]
    k = surviving[j]
    vkey = (str(level), GL, str(O(pi)), k)
    if vkey in Hecke_vector_cache:
        total = Hecke_vector_cache[vkey]
    else:
        reps, index_map = Hecke_transversal(pi)
        total = Hecke_operator_vector(All_Mat_Gens(GL)[k], pi, GL, reduction, reps, index_map)
        Hecke_vector_cache[vkey] = total
    v = Reduce_vector(total, taken, taken_columns)
    return vector(Reduction_ring(reduction[8]), [v[i] for i in surviving])

def Left_eigenvector(pis, system, L, GL=False, reduction=None):
    #The left eigenvector for the same eigenvalue system.  With w on the left,
    #a_rho = (w T_rho)_j / w_j needs only COLUMN j of T_rho, so harvesting an
    #eigenvalue costs one generator's worth of work rather than all of them.
    reduction = Current_reduction(GL, reduction, None)
    r = len(reduction[3])
    V = VectorSpace(L, r)
    for i in range(len(pis)):
        M = Hecke_matrix(pis[i], GL, reduction).change_ring(L) - L(system[i])*identity_matrix(L, r)
        V = V.intersection(M.kernel())          #Sage's .kernel() is the LEFT kernel
    assert V.dimension() >= 1, 'no left eigenvector for this system'
    return V.basis()[0], V.dimension()

def Harvest_eigenvalues(entry, pis, new_pis, GL=False, reduction=None, fast=True):
    #Given a certified 1-dimensional eigenline, return {rho: a_rho} for as many
    #new primes as wanted.  fast=True uses the left eigenvector and a single
    #column per prime; fast=False builds the whole operator and applies it to the
    #right eigenvector, which is the honest but r times more expensive route.
    assert Certified_line(entry), 'eigenspace is not 1-dimensional -- decompose with more primes first'
    reduction = Current_reduction(GL, reduction, None)
    L, system, v = entry[2], [entry[3][str(p)] for p in pis], entry[4]
    out = {}
    if fast:
        #The left eigenvector is only a computational device for reading off
        #eigenvalues one column at a time; the eigenvector of interest, returned
        #by Hecke_eigensystems, is the RIGHT one, v.  The condition that makes
        #this valid is that the left eigenspace for this system is a LINE: it is
        #then preserved by every commuting T_rho, so w T_rho = a_rho w with the
        #same a_rho that v sees.  Note w*v = 0 is NOT a problem -- in a Jordan
        #block the left and right eigenvectors for the same eigenvalue are always
        #orthogonal, so testing the pairing would reject the non-semisimple case
        #that harvests perfectly well.
        w, wdim = Left_eigenvector(pis, system, L, GL, reduction)
        assert wdim == 1, 'the left eigenspace is %s-dimensional, so it does not pin down the system; use fast=False' % wdim
        j = [i for i in range(len(w)) if w[i] != 0][0]
        for rho in new_pis:
            col = Hecke_column(rho, j, GL, reduction).change_ring(L)
            out[str(rho)] = (w*col)/w[j]
    else:
        k = [i for i in range(len(v)) if v[i] != 0][0]
        for rho in new_pis:
            T = Hecke_matrix(rho, GL, reduction).change_ring(L)
            image = T*v
            a = image[k]/v[k]
            assert image == a*v, 'not an eigenvector -- the certificate was violated'
            out[str(rho)] = a
    return out

def Decompose_until_certified(GL=False, reduction=None, coefficients=None, bound=30, max_primes=8):
    #Add primes one at a time and stop as soon as every piece is a certified
    #1-dimensional eigenline.  The decomposition is REFINED in place rather than
    #rebuilt from scratch on each new prime.
    reduction = Current_reduction(GL, reduction, coefficients)
    K = Reduction_ring(reduction[8])
    #The primes dividing the level must ALWAYS be candidates: Good_primes only
    #factors rationals below `bound`, so a level of norm 43 would never appear
    #with the default bound of 30 -- and without it the trace class never
    #separates, so nothing certifies.
    candidates = Good_primes(bound, True)
    for (pi_, e_) in factor(O(prime_base)):
        if not any(O(pi_) == O(q) for q in candidates):
            candidates.append(O(pi_))
    def priority(pi):
        #always take the primes DIVIDING THE LEVEL first -- U_pi at the level is
        #what separates a class from its Galois conjugate, and without it the
        #trace class stays merged -- then the primes above 2, whose behaviour is
        #exceptional (this is where the [D_u] eigenvector identity fails for
        #d = -1, -2, -7).  Everything else follows by norm, since N(pi) is the cost.
        divides_level = (O.fractional_ideal(level) + O.fractional_ideal(pi) == O.fractional_ideal(pi))
        even = (ZZ(O(pi).absolute_norm()) % 2 == 0)
        return (0 if divides_level else (1 if even else 2), ZZ(O(pi).absolute_norm()))
    candidates = sorted(candidates, key=priority)
    used = []
    pieces = [VectorSpace(K, len(reduction[3]))]
    systems = None
    for pi in candidates:
        used.append(pi)
        pieces = Refine_decomposition(pieces, pi, GL, reduction)
        systems = [[W.dimension()] + Eigenvalue_system(used, W, GL, reduction) for W in pieces]
        systems = [[e[0], e[4], e[1], dict(zip([str(p) for p in used], e[2])), e[3]] for e in systems]
        if all(Certified_line(e) for e in systems):
            return used, systems
        if len(used) >= max_primes:
            break
    return used, systems

def Level_is_Galois_stable():
    #Complex conjugation maps Gamma_0(p) to Gamma_0(conj p), so it is an
    #automorphism of Gamma_0(level) only when the level ideal is fixed by the
    #Galois action -- i.e. exactly when the prime base is RAMIFIED or INERT.
    I = O.fractional_ideal(level)
    return I == O.fractional_ideal(Conjugate_elt(level))

def Conjugate_matrix(M):
    #entrywise complex conjugation
    return Matrix([[Conjugate_elt(M[0][0]), Conjugate_elt(M[0][1])],
                   [Conjugate_elt(M[1][0]), Conjugate_elt(M[1][1])]])

def Fricke_element():
    return Matrix([[0,-1],[level,0]])

def Involution_matrix(which, GL=False, reduction=None, coefficients=None):
    #The matrix of an involution of Gamma_0(level) on the abelianisation, in the
    #basis of surviving generators.  which = 'fricke' is conjugation by
    #W = [[0,-1],[level,0]]: W M W^-1 = [[d, -c/level],[-level*b, a]], integral
    #because level divides c, with lower-left again divisible by level, so it
    #does preserve Gamma_0(level).  W^2 = -level*I is central, hence an involution
    #on the group.  which = 'conjugation' is entrywise complex conjugation, only
    #defined when the level is Galois stable.
    #These need NOT commute with the Hecke algebra, so they are analysed after the
    #eigenspaces are known rather than folded into the decomposition.
    reduction = Current_reduction(GL, reduction, coefficients)
    taken, taken_columns, surviving, n = reduction[1], reduction[2], reduction[3], reduction[4]
    generators = All_Mat_Gens(GL)
    to_vector = GElement_to_vector if GL else Element_to_vector
    if which == 'fricke':
        W = Fricke_element()
        Winv = W.inverse()
        act = lambda g: W*g*Winv
    elif which == 'conjugation':
        assert Level_is_Galois_stable(), 'complex conjugation needs a Galois-stable level: the prime base must be ramified or inert, not split'
        act = Conjugate_matrix
    else:
        raise ValueError("which must be 'fricke' or 'conjugation'")
    columns = []
    for k in surviving:
        image = act(generators[k])
        assert all(entry in O for entry in image.list()), 'the involution took a generator out of the ring'
        assert O(image[1][0]) in O.fractional_ideal(level), 'the involution left Gamma_0(level)'
        v = Reduce_vector(to_vector(image, n), taken, taken_columns)
        columns.append([v[j] for j in surviving])
    return Matrix(Reduction_ring(reduction[8]), columns).transpose()

def Fricke_matrix(GL=False, reduction=None, coefficients=None):
    return Involution_matrix('fricke', GL, reduction, coefficients)

def Conjugation_matrix(GL=False, reduction=None, coefficients=None):
    return Involution_matrix('conjugation', GL, reduction, coefficients)

def Hecke_system_of(w, pis, L, GL=False, reduction=None):
    #The Hecke eigenvalue system of a vector, or None if it is not a Hecke
    #eigenvector at all.
    reduction = Current_reduction(GL, reduction, None)
    k = [i for i in range(len(w)) if w[i] != 0]
    if not k:
        return None
    k = k[0]
    system = []
    for pi in pis:
        T = Hecke_matrix(pi, GL, reduction).change_ring(L)
        image = T*w
        a = image[k]/w[k]
        if image != a*w:
            return None
        system.append(a)
    return system

def Field_automorphisms(L):
    #automorphisms of the residue field: Galois maps for a number field, Frobenius
    #powers for a finite field, nothing to do for the prime field
    try:
        if L.degree() == 1:
            return []
    except Exception:
        return []
    try:
        return list(L.automorphisms())
    except Exception:
        pass
    try:
        return [L.frobenius_endomorphism(i) for i in range(L.degree())]
    except Exception:
        return []

def Analyse_involution(entry, systems, pis, which, GL=False, reduction=None):
    #Apply the involution to a Hecke eigenvector and report, in words, what
    #happened.  The primes are split into those COPRIME to the level and those
    #dividing it, because the two behave differently: the Fricke involution
    #commutes with T_pi away from the level but NOT with U_pi at it, so it can
    #preserve the away-from-level eigenvalue system while moving the vector
    #inside that system's eigenspace.  Reporting a single verdict over all primes
    #would just say "not an eigenvector" and hide that structure.
    reduction = Current_reduction(GL, reduction, None)
    L, v = entry[2], entry[4]
    F = Involution_matrix(which, GL, reduction).change_ring(L)
    w = F*v
    name = 'Fricke' if which == 'fricke' else 'complex conjugation'
    if w == 0:
        return '%s kills this eigenvector' % name, None, None
    k = [i for i in range(len(v)) if v[i] != 0][0]
    lam = w[k]/v[k]
    if w == lam*v:
        return 'IS an eigenvector under %s, with eigenvalue %s' % (name, lam), lam, None
    away  = [pi for pi in pis if O.fractional_ideal(level) + O.fractional_ideal(pi) != O.fractional_ideal(pi)]
    atlvl = [pi for pi in pis if pi not in away]
    notes = ['is NOT an eigenvector under %s' % name]
    #does the image keep the same eigenvalue system away from the level?
    #v's own eigenvalues, computed directly: entry[3] only holds the primes the
    #decomposition happened to use, which need not cover the analysis set
    sys_v = Hecke_system_of(v, away, L, GL, reduction)
    sys_away = Hecke_system_of(w, away, L, GL, reduction)
    if sys_away is not None and sys_v is not None and sys_away == sys_v:
        notes.append('but its image has the SAME T_pi system for every prime away from the level, so it stays inside that eigenspace')
    elif sys_away is not None:
        matched = None
        for j in range(len(systems)):
            try:
                if all(L(systems[j][3][str(away[i])]) == sys_away[i] for i in range(len(away))):
                    matched = j
                    break
            except Exception:
                continue
        if matched is not None:
            notes.append('its image lies in the eigenspace of piece %s (matching away from the level)' % matched)
        else:
            notes.append('its image has the system %s away from the level, which is not among the listed pieces'
                         % ', '.join(['a(%s)=%s' % (away[i], sys_away[i]) for i in range(len(away))]))
    else:
        notes.append('and its image is not a Hecke eigenvector even away from the level')
    #and what happened at the level itself
    for pi in atlvl:
        U = Hecke_matrix(pi, GL, reduction).change_ring(L)
        image = U*w
        j = [i for i in range(len(w)) if w[i] != 0][0]
        mu = image[j]/w[j]
        if image == mu*w:
            before = Hecke_system_of(v, [pi], L, GL, reduction)
            notes.append('at the level prime %s the eigenvalue moved from %s to %s'
                         % (pi, before[0] if before else '?', mu))
        else:
            notes.append('at the level prime %s the image is no longer a U_pi eigenvector' % pi)
    #is it the Galois conjugate of the same eigenvector?
    for sigma in Field_automorphisms(L):
        sv = vector(L, [sigma(c) for c in v])
        if sv == v:
            continue
        j = [i for i in range(len(sv)) if sv[i] != 0][0]
        mu = w[j]/sv[j]
        if w == mu*sv:
            notes.append('in fact the image is the GALOIS CONJUGATE eigenvector (%s -> %s), scaled by %s'
                         % (L.gen(), sigma(L.gen()), mu))
            break
    return '; '.join(notes), None, None

def Prime_base_above(rational_prime, which=0):
    #A generator of a prime of O above the given RATIONAL prime.  Split and
    #ramified primes give a proper element of O (not the rational prime itself);
    #`which` picks between the two primes above a split prime.  In the INERT case
    #the prime base must be handed back as a rational integer, because the inert
    #branch of Quotient_by_Level indexes residues by it.
    factors = list(factor(O(rational_prime)))
    if len(factors) == 1 and factors[0][1] == 1 and ZZ(factors[0][0].absolute_norm()) == rational_prime^2:
        return ZZ(rational_prime), 'inert'
    #factor() of an element of O returns ELEMENTS, not ideals
    element, e = factors[which % len(factors)]
    kind = 'ramified' if e > 1 else 'split'
    return O(element), kind

def Setup_group(alpha_value, rational_prime, exponent_value=1, which=0):
    #Build every global the rest of the notebook depends on, for one Bianchi group
    #and one level.  This is exactly what the setup cell does by hand.
    global alpha, half_unit_ord, x, O, f, prime_base, exponent, level
    global gen1, gen2, S, Fixed_Quotient, Real_Gens, Real_Rels, Real_Mat_Gens
    global Inv_gens, ZeroDiv_gens, zero_div
    alpha = alpha_value
    half_unit_ord = 2 if alpha == -1 else (3 if alpha == -3 else 1)
    x = fieldconstructor(alpha)
    O = (x.parent()).ring_of_integers()
    f = x.coordinates_in_terms_of_powers()
    prime_base, kind = Prime_base_above(rational_prime, which)
    exponent = exponent_value
    level = prime_base ^ exponent
    gen1 = matrix([[1,0],[1,1]])
    gen2 = matrix([[1,0],[x,1]])
    S = matrix([[0,-1],[1,0]])
    Fixed_Quotient = Quotient_by_Level()
    Real_Gens = special_generators()
    Real_Rels = Reidemeister_Schreier_Relations()
    Real_Mat_Gens = Matrix_Generators()
    Inv_gens = Seperate_Inv_generators(Real_Gens)
    ZeroDiv_gens = Seperate_PriPow_generators(Real_Gens)
    zero_div = ZeroDivsQuotient()
    Clear_Hecke_cache()
    return kind

def Hecke_algebra(GL=False, Projective=False, alpha=None, prime=2, exponent=1,
                  coefficients=0, which=0, bound=30, harvest=0, verbose=True):
    #The whole pipeline for one group and one level.
    #  GL           False -> SL2, True -> GL2
    #  Projective   False -> SL/GL,  True -> PSL/PGL
    #  alpha        the Bianchi field; None keeps the current global
    #  prime        a RATIONAL prime; prime_base is a prime of O above it, which
    #               is not the rational prime itself unless it is inert
    #  exponent     level = prime_base^exponent
    #  coefficients 0 for the free part, a prime p for the p-torsion part
    #  which        which prime above `prime` to take, when it splits
    #  harvest      extra eigenvalues to collect per piece, beyond the ones the
    #               decomposition needed (0 = none)
    #Returns a dictionary; with verbose=True it also prints a readable report.
    global Fixed_Reduction
    if alpha is None:
        alpha = globals()['alpha']
    kind = Setup_group(alpha, prime, exponent, which)
    assert coefficients == 0 or is_prime(coefficients), 'coefficients must be 0 (free part) or a rational prime (torsion part)'
    integral = Abelianisation_reduction(Projective, GL, -1)
    divisors = Elementary_divisors(Projective, GL, integral)
    tprimes = Torsion_primes(Projective, GL, integral)
    free_rank = len(Abelianisation_reduction(Projective, GL, 0)[3])
    Fixed_Reduction = Abelianisation_reduction(Projective, GL, coefficients)
    dim = len(Fixed_Reduction[3])
    group = ('P' if Projective else '') + ('GL' if GL else 'SL')
    out = {'group': group, 'alpha': alpha, 'prime': prime, 'prime_base': prime_base,
           'kind': kind, 'level': level, 'norm': ZZ(O(level).absolute_norm()),
           'elementary_divisors': divisors, 'torsion_primes': tprimes,
           'phi_primes': Phi_primes(), 'free_rank': free_rank,
           'coefficients': coefficients, 'dimension': dim,
           'reduction': Fixed_Reduction, 'pieces': []}
    if verbose:
        print('%s_0(%s) inside %s_2(O_%s)   level = %s^%s = %s of norm %s (%s)'
              % (('P' if Projective else '') + ('GL' if GL else 'SL'), level, 'GL' if GL else 'SL',
                 alpha, prime_base, exponent, level, out['norm'], kind))
        print('  abelianisation : Z^%s + %s      torsion primes %s   (phi predicts %s)'
              % (free_rank, divisors, tprimes, out['phi_primes']))
        print('  working %s : dimension %s'
              % ('over QQ (free part)' if coefficients == 0 else 'over GF(%s) (%s-torsion)' % (coefficients, coefficients), dim))
        if coefficients != 0 and coefficients not in tprimes:
            print('  NOTE: %s does not divide the torsion, so this is just the free part reduced mod %s' % (coefficients, coefficients))
    if dim == 0:
        if verbose:
            print('  nothing to decompose at these coefficients')
        return out
    used, systems = Decompose_until_certified(GL, Fixed_Reduction, coefficients, bound)
    out['decomposition_primes'] = used
    out['unit_hecke'] = Unit_Hecke_matrix(GL, Fixed_Reduction)
    out['fricke'] = Fricke_matrix(GL, Fixed_Reduction)
    out['galois_stable'] = Level_is_Galois_stable()
    out['conjugation'] = Conjugation_matrix(GL, Fixed_Reduction) if out['galois_stable'] else None
    analysis_primes = sorted(set(list(used) + Primes_by_norm(14, True, 14)),
                             key=lambda q: ZZ(O(q).absolute_norm()))
    #Always collect a few eigenvalues at small primes AWAY from the level, even
    #when the decomposition did not need them: the report is uninformative
    #without them, and the Eisenstein/trace classification is meaningless at the
    #level prime (where the trace class has eigenvalue pi itself, not Tr(pi)).
    base = [pi for pi in Primes_by_norm(14, False, 14) if not any(pi == q for q in used)]
    extra = base + [pi for pi in Primes_by_norm(bound, True)
                    if not any(pi == q for q in used) and not any(pi == q for q in base)][:harvest]
    if verbose:
        print('  decomposed with %s into %s piece(s)' % ([str(p) for p in used], len(systems)))
        if not out['galois_stable']:
            print('  complex conjugation is not available: the level is split, so it is not Galois stable')
    for idx in range(len(systems)):
        e = systems[idx]
        L = e[2]
        record = {'generalised_dim': e[0], 'eigenspace_dim': e[1], 'field': L,
                  'eigenvalues': dict(e[3]), 'eigenvector': e[4],
                  'certified': Certified_line(e)}
        if extra and Certified_line(e):
            record['eigenvalues'].update(Harvest_eigenvalues(e, used, extra, GL, Fixed_Reduction, True))
        #classify on the primes AWAY from the level only: at the level the trace
        #class has eigenvalue pi itself, not Tr(pi), so including it would label
        #every trace class 'other'
        known = record['eigenvalues']
        away = [pi for pi in list(used) + extra
                if str(pi) in known and O.fractional_ideal(level) + O.fractional_ideal(pi) != O.fractional_ideal(pi)]
        allE = bool(away) and all(known[str(pi)] == L(ZZ(pi.absolute_norm())+1) for pi in away)
        allT = bool(away) and all(known[str(pi)] == L(ZZ(pi.trace())) for pi in away)
        record['type'] = 'Eisenstein' if allE else ('trace' if allT else 'other')
        record['fricke'] = Analyse_involution(e, systems, analysis_primes, 'fricke', GL, Fixed_Reduction)[0]
        if out['galois_stable']:
            record['conjugation'] = Analyse_involution(e, systems, analysis_primes, 'conjugation', GL, Fixed_Reduction)[0]
        out['pieces'].append(record)
        if verbose:
            semis = 'semisimple' if e[0] == (L.degree() if hasattr(L,'degree') else 1)*e[1] else 'JORDAN BLOCK'
            cert = '' if record['certified'] else '   NOT CERTIFIED (eigenspace dim > 1: add more primes, or raise bound/max_primes)'
            print('  piece %s : generalised dim %s, eigenspace dim %s, %s, over %s   [%s]%s'
                  % (idx, e[0], e[1], semis, L, record['type'], cert))
            keys = sorted(record['eigenvalues'].keys())
            print('     %s' % '  '.join(['a(%s)=%s' % (k, record['eigenvalues'][k]) for k in keys]))
            print('     Fricke      : %s' % record['fricke'])
            if out['galois_stable']:
                print('     conjugation : %s' % record['conjugation'])
    return out

def Structure_string(free_rank, divisors):
    #'Z^2 + Z/2 + Z/10' -- the invariant factors as Sage's smith_form orders them
    parts = []
    if free_rank == 1:
        parts.append('Z')
    elif free_rank > 1:
        parts.append('Z^%s' % free_rank)
    parts += ['Z/%s' % d for d in divisors]
    return ' + '.join(parts) if parts else '0'

def Subgroup_index():
    #[SL_2(O) : Gamma_0(level)] = |P^1(O/level)| = N(p)^(e-1) * (N(p)+1)
    N = ZZ(O(level).absolute_norm())
    p_norm = ZZ(O(prime_base).absolute_norm())
    return (N//p_norm)*(p_norm + 1)

def Abelianisation(GL=False, Projective=False, alpha=None, prime=2, exponent=1,
                   which=0, verbose=True):
    #The abelianisation of Gamma_0(level), for one Bianchi group and one level.
    #Same inputs as Hecke_algebra minus the coefficient ring, since this is the
    #full integral answer: the free rank and the torsion together.
    #  GL           False -> SL2, True -> GL2
    #  Projective   False -> SL/GL,  True -> PSL/PGL
    #  alpha        the Bianchi field; None keeps the current global
    #  prime        a RATIONAL prime; prime_base is a prime of O above it
    #  exponent     level = prime_base^exponent
    #  which        which prime above `prime` to take, when it splits
    #Returns a dictionary and, with verbose=True, prints a readable report.
    global Fixed_Reduction
    if alpha is None:
        alpha = globals()['alpha']
    kind = Setup_group(alpha, prime, exponent, which)
    Fixed_Reduction = Abelianisation_reduction(Projective, GL, -1)
    reduced, surviving = Fixed_Reduction[0], Fixed_Reduction[3]
    divisors = Elementary_divisors(Projective, GL, Fixed_Reduction)
    #free rank straight from the reduced lattice, avoiding a second matrix build
    free_rank = len(surviving) - (reduced.rank() if reduced.nrows() else 0)
    tprimes = Torsion_primes(Projective, GL, Fixed_Reduction)
    group = ('P' if Projective else '') + ('GL' if GL else 'SL')
    out = {'group': group, 'alpha': alpha, 'prime': prime, 'prime_base': prime_base,
           'exponent': exponent, 'kind': kind, 'level': level,
           'norm': ZZ(O(level).absolute_norm()), 'index': Subgroup_index(),
           'free_rank': free_rank, 'elementary_divisors': divisors,
           'structure': Structure_string(free_rank, divisors),
           'torsion_primes': tprimes, 'phi_primes': Phi_primes(),
           'surviving_generators': surviving, 'reduction': Fixed_Reduction,
           'relevant_coefficients': [0] + tprimes}
    if verbose:
        print('%s_0(%s) inside %s_2(O_%s)   level = (%s)^%s of norm %s (%s)   index %s'
              % (group, level, 'GL' if GL else 'SL', alpha, prime_base, exponent,
                 out['norm'], kind, out['index']))
        print('  abelianisation : %s' % out['structure'])
        print('  free rank %s, elementary divisors %s, %s surviving generators'
              % (free_rank, divisors, len(surviving)))
        print('  torsion primes %s   (phi predicts %s)%s'
              % (tprimes, out['phi_primes'],
                 '   <== %s not predicted by phi' % [q for q in tprimes if q not in out['phi_primes']]
                 if [q for q in tprimes if q not in out['phi_primes']] else ''))
        print('  coefficient rings worth running through Hecke_algebra: %s' % out['relevant_coefficients'])
    return out


# ----------------------------------------------------------------------------
#  interactive front end
# ----------------------------------------------------------------------------

EUCLIDEAN = [-1, -2, -3, -7, -11]

SUPERSCRIPT = {'0':'⁰','1':'¹','2':'²','3':'³','4':'⁴',
               '5':'⁵','6':'⁶','7':'⁷','8':'⁸','9':'⁹'}

def superscript(n):
    return ''.join(SUPERSCRIPT.get(c, c) for c in str(n))

def root(d):
    return '√%s' % d

def generator_description(d):
    if d % 4 == 1:
        return '(1 + %s)/2' % root(d)
    return root(d)

def ask_field():
    print('')
    print('(1) Which Bianchi group?')
    print('    The five Euclidean cases are  d = -1, -2, -3, -7, -11')
    while True:
        answer = input('    d = ').strip()
        try:
            d = Integer(answer)
        except Exception:
            print('        "%s" is not an integer.  Try again.' % answer)
            continue
        if d in EUCLIDEAN:
            print('        working in O = Z[x] ⊂ Q(%s),  with  x = %s'
                  % (root(d), generator_description(d)))
            return d
        print('        %s is not one of the five Euclidean discriminants.  Try again.' % d)

def splitting_type(d, p):
    x = fieldconstructor(d)
    K = x.parent()
    factors = list(K.ideal(p).factor())
    if len(factors) == 2:
        return 'split', ZZ(factors[0][0].absolute_norm())
    if factors[0][1] > 1:
        return 'ramified', ZZ(factors[0][0].absolute_norm())
    return 'inert', ZZ(factors[0][0].absolute_norm())

def ask_prime(d):
    print('')
    print('(2) Which rational prime divides the norm of the level?')
    while True:
        answer = input('    p = ').strip()
        try:
            p = Integer(answer)
        except Exception:
            print('        "%s" is not an integer.  Try again.' % answer)
            continue
        if p > 1 and is_prime(p):
            kind, norm = splitting_type(d, p)
            print('        %s is %s in Q(%s)' % (p, kind, root(d)))
            return p, norm, kind
        print('        %s is not a rational prime.  Try again.' % p)

def ask_exponent(prime_norm):
    print('')
    print('(3) What is the exponent of the level?   (level = πᵉ, with π the prime above p)')
    while True:
        answer = input('    e = ').strip()
        try:
            e = Integer(answer)
        except Exception:
            print('        "%s" is not an integer.  Try again.' % answer)
            continue
        if e >= 1:
            index = (prime_norm^(e-1))*(prime_norm + 1)
            print('        Γ₀(level) has index %s in the full Bianchi group  (level norm %s)'
                  % (index, prime_norm^e))
            return e
        print('        The exponent must be at least 1.  Try again.')

def ask_group():
    print('')
    print('(4) Which group?   SL, PSL, GL, PGL')
    table = {'SL': (False, False), 'PSL': (True, False),
             'GL': (False, True),  'PGL': (True, True)}
    while True:
        answer = input('    group = ').strip().upper()
        if answer in table:
            return table[answer], answer
        print('        "%s" is not one of SL, PSL, GL, PGL.  Try again.' % answer)

def ask_coefficients():
    print('')
    print('(5) Which coefficients for the Hecke action?')
    print('    free   the free part; eigenvalues in number fields, exact (not floating point)')
    print('    p      a rational prime, for the mod p action; eigenvalues in GF(p^k)')
    while True:
        answer = input('    coefficients = ').strip().lower()
        if answer in ('free', 'f', '0', 'q', 'qq', 'c', 'cc'):
            return 0
        try:
            q = Integer(answer)
        except Exception:
            print('        type "free", or a rational prime.  Try again.')
            continue
        if q > 1 and is_prime(q):
            return q
        print('        %s is neither "free" nor a rational prime.  Try again.' % answer)

def ask_how_many():
    print('')
    print('(6) How many Hecke operators should be computed?')
    print('    In order to determine Hecke eigenvalue systems efficiently, Hecke')
    print('    operators indexed by even primes and primes dividing the level are')
    print('    always computed.  When the group is SL or PSL the "unit Hecke operator"')
    print('    is also included.')
    print('    Save for possibly the one at the level, the operators are computed in')
    print('    ascending order of the norm of the prime indexing them.')
    while True:
        answer = input('    how many = ').strip()
        try:
            n = Integer(answer)
        except Exception:
            print('        "%s" is not an integer.  Try again.' % answer)
            continue
        if n >= 1:
            return n
        print('        ask for at least 1.  Try again.')

def ask_extra_primes():
    print('')
    print('(7) Any particular rational primes you want Hecke operators for?')
    print('    Separate them by spaces or commas; leave blank for none.')
    print('    (Every prime of O above each one you name will be used.)')
    while True:
        answer = input('    primes = ').strip()
        if answer == '':
            return []
        tokens = [t for t in answer.replace(',', ' ').split() if t]
        out = []
        bad = None
        for t in tokens:
            try:
                q = Integer(t)
            except Exception:
                bad = t; break
            if not (q > 1 and is_prime(q)):
                bad = t; break
            out.append(q)
        if bad is None:
            return out
        print('        "%s" is not a rational prime.  Try again (or leave blank).' % bad)

def hecke_operator_primes(count, extra_rationals):
    #Ordering: the primes dividing the LEVEL first, then the primes above 2, then
    #everything else by ascending norm.  The first two groups are always included
    #whatever `count` is, because they are what actually separates the eigenspaces.
    candidates = Good_primes(80, True)
    for (pi_, e_) in factor(O(prime_base)):
        if not any(O(pi_) == O(q) for q in candidates):
            candidates.append(O(pi_))
    def divides_level(pi):
        return (O.fractional_ideal(level) + O.fractional_ideal(pi)
                == O.fractional_ideal(pi))
    def priority(pi):
        n = ZZ(O(pi).absolute_norm())
        return (0 if divides_level(pi) else (1 if n % 2 == 0 else 2), n)
    ordered = sorted(candidates, key=priority)
    mandatory = [pi for pi in ordered
                 if divides_level(pi) or ZZ(O(pi).absolute_norm()) % 2 == 0]
    chosen = list(mandatory)
    for pi in ordered:
        if len(chosen) >= count:
            break
        if not any(pi == q for q in chosen):
            chosen.append(pi)
    for q in extra_rationals:
        for (pi_, e_) in factor(O(q)):
            if not any(O(pi_) == O(c) for c in chosen):
                chosen.append(O(pi_))
    return sorted(chosen, key=priority), len(mandatory)

def field_definition(L):
    #how the generator of the residue field is defined
    try:
        if L.degree() == 1:
            return None
    except Exception:
        return None
    try:
        return '%s is a root of  %s' % (L.gen(), L.modulus())          # finite field
    except Exception:
        pass
    try:
        return '%s is a root of  %s' % (L.gen(), L.defining_polynomial())  # number field
    except Exception:
        return None

def piece_containing(w, pieces, L):
    #Which piece's span contains the vector w?  Complex conjugation sends the
    #eigenspace of the system (a_pi) to that of (a_conj(pi)), so it permutes the
    #pieces; when it does not fix one, the image lands inside another, and this
    #finds which.  The pieces live over the base field, so they are base-changed
    #to L before testing membership.
    for j in range(len(pieces)):
        W = pieces[j]
        if W.dimension() == 0:
            continue
        span = W.basis_matrix().change_ring(L).row_space()
        if w in span:
            return j
    return None

def operator_eigenvalue(M, v, L):
    #the scalar by which M acts on the eigenline through v, or None if it moves it
    w = M.change_ring(L)*v
    nz = [i for i in range(len(v)) if v[i] != 0]
    if not nz:
        return None
    k = nz[0]
    lam = w[k]/v[k]
    if w == lam*v:
        return lam
    return None

def classify(system, test, L):
    #Eisenstein and trace are judged AWAY FROM THE LEVEL: at the level the trace
    #class has eigenvalue pi itself, not Tr(pi), so including it would mislabel.
    away = [i for i in range(len(test))
            if O.fractional_ideal(level) + O.fractional_ideal(test[i])
               != O.fractional_ideal(test[i])]
    if not away:
        return 'unclassified (only the level operator was computed)'
    if all(system[i] == L(ZZ(O(test[i]).absolute_norm())+1) for i in away):
        return 'EISENSTEIN   a(π) = N(π)+1 for every π NOT dividing the level'
    if all(system[i] == L(ZZ(O(test[i]).trace())) for i in away):
        return 'trace        a(π) = Tr(π) for every π NOT dividing the level'
    return 'neither Eisenstein nor trace'

def main():
    print('=' * 76)
    print(' Hecke algebra on H₁ of Γ₀ congruence subgroups of Bianchi groups')
    print(' version %s   core sha256:%s' % (SCRIPT_VERSION, CORE_FINGERPRINT))
    print('=' * 76)
    d = ask_field()
    p, prime_norm, kind = ask_prime(d)
    e = ask_exponent(prime_norm)
    (projective, general), label = ask_group()
    coefficients = ask_coefficients()
    count = ask_how_many()
    extra = ask_extra_primes()

    print('')
    print('    building the group and the relation matrix ...')
    Setup_group(d, p, e, 0)
    #ONE relation-matrix build, reduced directly to the requested coefficients.
    #(The library's Hecke_algebra also computes the integral and rational
    #reductions for its report, rebuilding the matrix each time; the build is the
    #dominant cost, so doing it once here is about three times faster.)
    global Fixed_Reduction
    Fixed_Reduction = Abelianisation_reduction(projective, general, coefficients)
    dim = len(Fixed_Reduction[3])
    base = 'Q' if coefficients == 0 else 'F_%s' % coefficients
    ring = 'the free part over Q' if coefficients == 0 else 'mod %s' % coefficients

    print('')
    print('-' * 76)
    print(' Γ₀(%s) ⊂ %s₂(O),   O = Z[x] with x = %s'
          % (level, label, generator_description(d)))
    print(' level = (%s)%s of norm %s,  %s,  index %s'
          % (prime_base, superscript(e), ZZ(O(level).absolute_norm()), kind, Subgroup_index()))
    print(' coefficients: %s        dimension %s' % (ring, dim))
    if dim == 0:
        print('')
        print('   that space is zero -- no Hecke action to describe.')
        #Nothing else is computed on this path, so it is worth one more relation
        #matrix build to say which coefficients WOULD give something, instead of
        #leaving the user to guess and re-run.
        integral = Abelianisation_reduction(projective, general, -1)
        divisors = Elementary_divisors(projective, general, integral)
        torsion = Torsion_primes(projective, general, integral)
        free_rank = len(integral[3]) - (integral[0].rank() if integral[0].nrows() else 0)
        print('')
        print('   the abelianisation is %s' % Structure_string(free_rank, divisors))
        print('   free rank %s,  torsion primes %s' % (free_rank, torsion))
        options = ([] if free_rank == 0 else ['free']) + [str(q) for q in torsion]
        if coefficients == 0:
            print('   there is no free part, so "free" sees nothing.')
        else:
            print('   %s divides none of the invariant factors, so it sees nothing.' % coefficients)
        if options:
            print('   coefficients that do see something here: %s' % ', '.join(options))
        else:
            print('   the abelianisation is trivial: no choice of coefficients sees anything.')
        print('-' * 76)
        return None

    test, n_mandatory = hecke_operator_primes(count, extra)
    #The unit Hecke operator is conjugation by diag(1,u): u = x for d = -1,-3 and
    #u = -1 otherwise, so it is written T_x or T_-1 to match T_pi.  It is listed
    #FIRST, with N = 1, since the unit generates the trivial ideal.  (In GL and
    #PGL that conjugation is inner, so the operator is the identity there and
    #there is nothing to report.)
    unit_matrix_op = None
    unit_label = None
    if not general:
        unit_label = '%s' % Unit_generator()
        unit_matrix_op = Unit_Hecke_matrix(general, Fixed_Reduction)
    total = len(test) + (1 if unit_matrix_op is not None else 0)
    #"even prime(s)": 2 splits only for d = -7, so in four of the five fields
    #there is exactly one prime above 2.  When the level itself sits above 2 the
    #two categories overlap -- the count is a union so it does not double count,
    #but the wording must not claim they are separate.  (For d = -7 and p = 2 the
    #level is one of TWO even primes, hence "include" rather than "is".)
    level_is_even = (ZZ(O(prime_base).absolute_norm()) % 2 == 0)
    if level_is_even:
        mandatory_text = ('the unit and the even prime(s), which include the level'
                          if unit_matrix_op is not None
                          else 'the even prime(s), which include the level')
    else:
        mandatory_text = ('the unit, the level and the even prime(s)'
                          if unit_matrix_op is not None
                          else 'the level and the even prime(s)')
    print('')
    print(' computing %s Hecke operators (%s mandatory: %s):'
          % (total, n_mandatory + (1 if unit_matrix_op is not None else 0), mandatory_text))
    if unit_matrix_op is not None:
        print('     T_%-11s  N = %-5s (unit operator)' % (unit_label, 1))
    for q in test:
        tag = []
        if O.fractional_ideal(level) + O.fractional_ideal(q) == O.fractional_ideal(q):
            tag.append('at the level')
        if ZZ(O(q).absolute_norm()) % 2 == 0:
            tag.append('even')
        print('     T_%-11s  N = %-5s %s' % (q, ZZ(O(q).absolute_norm()),
                                             ('(%s)' % ', '.join(tag)) if tag else ''))
    #The Fricke involution is conjugation by W = [[0,-1],[level,0]].  It is
    #defined at every level and is cheap (one decomposition per generator, not
    #per generator per coset), so it is always computed.
    print('     Fricke involution')
    fricke_op = Fricke_matrix(general, Fixed_Reduction)
    #Complex conjugation is an automorphism of Gamma_0(level) exactly when the
    #level ideal is Galois stable, i.e. when the prime is inert or ramified; for
    #a split prime it maps Gamma_0(pi) to Gamma_0(pibar) and there is nothing to
    #report.  No need to ask -- just compute it where it exists.
    conj_op = None
    if kind in ('inert', 'ramified'):
        print('     complex conjugation')
        conj_op = Conjugation_matrix(general, Fixed_Reduction)

    print('')
    print(' characteristic polynomials, factored over %s:' % base)
    if unit_matrix_op is not None:
        print('     T_%-11s  %s' % (unit_label, unit_matrix_op.charpoly().factor()))
    for q in test:
        T = Hecke_matrix(q, general, Fixed_Reduction)
        print('     T_%-11s  %s' % (q, T.charpoly().factor()))
    print('     Fricke        %s' % fricke_op.charpoly().factor())
    if conj_op is not None:
        print('     conjugation   %s' % conj_op.charpoly().factor())

    print('')
    print('    decomposing ...')
    pieces = Hecke_decomposition(test, general, Fixed_Reduction)
    resolved = []
    for W in pieces:
        try:
            resolved.append(Eigenvalue_system(test, W, general, Fixed_Reduction))
        except Exception as err:
            resolved.append(None)
    print('')
    print(' %s Galois orbit(s) of eigenvalue systems:' % len(pieces))
    for idx in range(len(pieces)):
        W = pieces[idx]
        if resolved[idx] is None:
            print('')
            print('   piece %s:  dimension %s -- UNRESOLVED' % (idx+1, W.dimension()))
            print('              ask for more Hecke operators to split it.')
            continue
        L, system, eigvec, edim = resolved[idx]
        deg = L.degree() if hasattr(L, 'degree') else 1
        print('')
        #the generalised eigenspace is the primary component; over the residue
        #field its dimension is the piece's dimension divided by the orbit size.
        #It equals the eigenspace dimension exactly when the action is semisimple.
        gen_dim = W.dimension()//deg
        print('   piece %s:  dimension %s,  eigenspace dimension %s,  generalised eigenspace dimension %s,'
              % (idx+1, W.dimension(), edim, gen_dim))
        print('              residue field %s' % L)
        if deg > 1:
            print('              contains %s CONJUGATE eigenclasses (one Galois orbit over %s);'
                  % (deg, base))
            print('              the system below is one of them, the others are its')
            print('              images under the Galois action of %s over %s' % (L, base))
        else:
            print('              a single eigenclass, defined over %s' % base)
        definition = field_definition(L)
        if definition is not None:
            print('              where %s' % definition)
        print('              %s' % classify(system, test, L))
        if edim > 1:
            print('              NOT CERTIFIED: eigenspace dimension > 1, so a further')
            print('              operator could split this piece.')
        lam = operator_eigenvalue(fricke_op, eigvec, L)
        if lam is not None:
            print('              Fricke involution PRESERVES this eigenclass, acting by %s' % lam)
        else:
            print('              Fricke involution does not preserve this eigenclass')
        if conj_op is not None:
            lam = operator_eigenvalue(conj_op, eigvec, L)
            if lam is not None:
                print('              complex conjugation PRESERVES this eigenclass, acting by %s' % lam)
            else:
                w = conj_op.change_ring(L)*eigvec
                target = piece_containing(w, pieces, L)
                if target is None:
                    print('              complex conjugation moves this eigenclass, but its image')
                    print('              lies in no single piece -- add more Hecke operators')
                elif target == idx:
                    print('              complex conjugation preserves the PIECE but not the')
                    print('              eigenline: it sends this class to a Galois conjugate')
                    print('              of itself inside piece %s' % (idx+1))
                else:
                    print('              complex conjugation sends this eigenclass to piece %s:'
                          % (target+1))
                    print('              they form a conjugate pair')
        if unit_matrix_op is not None:
            lam = operator_eigenvalue(unit_matrix_op, eigvec, L)
            print('                 a(%-11s) = %-24s    <- unit operator'
                  % (unit_label, lam if lam is not None else 'no scalar'))
        for k in range(len(test)):
            print('                 a(%-11s) = %s' % (test[k], system[k]))
    print('-' * 76)
    return pieces

# Sage's script runner sets __name__ to "sage.all" rather than "__main__", so the
# usual guard would never fire.  Call main() directly.
main()
