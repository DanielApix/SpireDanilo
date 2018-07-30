from multiprocessing import Pool

def index_in_alphabet(t, typ_alphabet_list):
    return typ_alphabet_list.index(t)


# ------------------------ CFL ---------------------------------------------------------------------
# CFL - Lyndon factorization - Duval's algorithm
def CFL(word):
    """
    CFL Duval's algorithm.
    """
    CFL_list = []
    k = 0
    while k < len(word):
        i = k + 1
        j = k + 2
        while True:
            if j == len(word) + 1 or word[j - 1] < word[i - 1]:
                while k < i:
                    # print(word[k:k + j - i])
                    CFL_list.append(word[k:k + j - i])
                    k = k + j - i
                break
            else:
                # if word[j-1] > word[i-1]:
                if word[j - 1] > word[i - 1]:
                    i = k + 1
                else:
                    i = i + 1
                j = j + 1

    return CFL_list


# CFL - Lyndon factorization - Duval's algorithm - on a specific algorithm
def CFL_for_alphabet(word, list_alphabet):
    """
    CFL Duval's algorithm.
    """
    CFL_list = []
    k = 0
    while k < len(word):
        i = k + 1
        j = k + 2
        while True:
            if j == len(word) + 1 or index_in_alphabet(word[j - 1], list_alphabet) < index_in_alphabet(word[i - 1],
                                                                                                       list_alphabet):
                while k < i:
                    # print(word[k:k + j - i])
                    CFL_list.append(word[k:k + j - i])
                    k = k + j - i
                break
            else:
                # if word[j-1] > word[i-1]:
                if index_in_alphabet(word[j - 1], list_alphabet) > index_in_alphabet(word[i - 1], list_alphabet):
                    i = k + 1
                else:
                    i = i + 1
                j = j + 1

    return CFL_list


# ---------------------------------------------------------------------------------------------



# ----------------------- ICFL ----------------------------------------------------------------------

# ICFL recursive (without using of compute_br)- Inverse Lyndon factorization
def ICFL_recursive(word):
    """In this version of ICFL, we don't execute compute_br - one only O(n) scanning of word"""
    #br_list = []
    icfl_list = []

    #compute_icfl_recursive(word, br_list, icfl_list)
    compute_icfl_recursive(word, icfl_list)

    return icfl_list


def ICFL_recursive_for_alphabet(word, list_alphabet):
    """In this version of ICFL, we don't execute compute_br - one only O(n) scanning of word"""
    #br_list = []
    icfl_list = []

    #compute_icfl_recursive_for_alphabet(word, br_list, icfl_list, list_alphabet)
    compute_icfl_recursive_for_alphabet(word, icfl_list, list_alphabet)

    return icfl_list


def compute_icfl_recursive(word, icfl_list):
    # At each step compute the current bre
    pre_pair = find_pre(word)
    current_bre_quad = find_bre(pre_pair)
    #br_list.append(current_bre_quad)

    if current_bre_quad[1] == '' and current_bre_quad[0].find('$') >= 0:
        w = current_bre_quad[0]
        icfl_list.insert(0, w[:len(w) - 1])
        return
    else:
        compute_icfl_recursive(current_bre_quad[1] + current_bre_quad[2], icfl_list)
        if len(icfl_list[0]) > current_bre_quad[3]:
            icfl_list.insert(0, current_bre_quad[0])
        else:
            icfl_list[0] = current_bre_quad[0] + icfl_list[0]
        return


# ICFL recursive (without using of compute_br)- Inverse Lyndon factorization - on a specific algorithm
def compute_icfl_recursive_for_alphabet(word, icfl_list, list_alphabet):
    # At each step compute the current bre
    pre_pair = find_pre_for_alphabet(word, list_alphabet)
    current_bre_quad = find_bre_for_alphabet(pre_pair, list_alphabet)
    #br_list.append(current_bre_quad)

    if current_bre_quad[1] == '' and current_bre_quad[0].find('$') >= 0:
        w = current_bre_quad[0]
        icfl_list.insert(0, w[:len(w) - 1])
        return
    else:
        compute_icfl_recursive_for_alphabet(current_bre_quad[1] + current_bre_quad[2], icfl_list,
                                            list_alphabet)
        if len(icfl_list[0]) > current_bre_quad[3]:
            icfl_list.insert(0, current_bre_quad[0])
        else:
            icfl_list[0] = current_bre_quad[0] + icfl_list[0]
        return


def find_pre(word):
    if len(word) == 1:
        return (word + "$", '')
    else:
        i = 0
        j = 1

        while j < len(word) and word[j] <= word[i]:
            if word[j] < word[i]:
                i = 0
            else:
                i = i + 1
            j = j + 1

        if j == len(word):
            return (word + "$", '')
        else:
            return (word[0:j + 1], word[j + 1:len(word)])


def find_pre_for_alphabet(word, list_alphabet):
    if len(word) == 1:
        return (word + "$", '')
    else:
        i = 0
        j = 1

        while j < len(word) and index_in_alphabet(word[j], list_alphabet) <= index_in_alphabet(word[i], list_alphabet):
            if index_in_alphabet(word[j], list_alphabet) < index_in_alphabet(word[i], list_alphabet):
                i = 0
            else:
                i = i + 1
            j = j + 1

        if j == len(word):
            return (word + "$", '')
        else:
            return (word[0:j + 1], word[j + 1:len(word)])


def find_bre(pre_pair):
    w = pre_pair[0]
    v = pre_pair[1]

    if v == '' and w.find('$') >= 0:
        # return (w[:len(w)-1], '', '', 0)
        return (w, '', '', 0)
    else:
        n = len(w) - 1

        f = border(w[:n])

        i = n
        last = f[i - 1]

        while i > 0:
            if w[f[i - 1]] < w[n]:
                last = f[i - 1]
            i = f[i - 1]

        return (w[:n - last], w[n - last:n + 1], v, last)


def find_bre_for_alphabet(pre_pair, list_alphabet):
    w = pre_pair[0]
    v = pre_pair[1]

    if v == '' and w.find('$') >= 0:
        # return (w[:len(w)-1], '', '', 0)
        return (w, '', '', 0)
    else:
        n = len(w) - 1

        f = border(w[:n])

        i = n
        last = f[i - 1]

        while i > 0:
            if index_in_alphabet(w[f[i - 1]], list_alphabet) < index_in_alphabet(w[n], list_alphabet):
                last = f[i - 1]
            i = f[i - 1]

        return (w[:n - last], w[n - last:n + 1], v, last)


def border(p):
    l = len(p)
    pi = [0]
    k = 0
    for i in range(1, l):
        while (k > 0 and p[k] != p[i]):
            k = pi[k - 1]
        if (p[k] == p[i]):
            pi.append(k + 1)
            k = k + 1
        else:
            pi.append(k)

    return pi



# ------------------------ CFL_icfl ---------------------------------------------------------------------
# CFL factorization - ICFL subdecomposition
def CFL_icfl(word, C):
    """
    CFL Duval's algorithm.
    """
    CFL_list = []
    k = 0
    while k < len(word):
        i = k + 1
        j = k + 2
        while True:
            if j == len(word) + 1 or word[j - 1] < word[i - 1]:
                while k < i:
                    # print(word[k:k + j - i])
                    w = word[k:k + j - i]
                    if len(w) <= C:
                        CFL_list.append(word[k:k + j - i])
                    else:
                        ICFL_list_recursive = ICFL_recursive(w)
                        p= Pool(5)

                        # Insert << to indicate the begin of the subdecomposition of w
                        CFL_list.append("<<")
                        for v in ICFL_list_recursive:
                            CFL_list.append(v)
                        # Insert >> to indicate the end of the subdecomposition of w
                        CFL_list.append(">>")

                    k = k + j - i
                break
            else:
                # if word[j-1] > word[i-1]:
                if word[j - 1] > word[i - 1]:
                    i = k + 1
                else:
                    i = i + 1
                j = j + 1

    return CFL_list


# CFL factorization - ICFL subdecomposition - on a specific algorithm
def CFL_icfl_for_alphabet(word, C, list_alphabet):
    """
    CFL Duval's algorithm.
    """
    CFL_list = []
    k = 0
    while k < len(word):
        i = k + 1
        j = k + 2
        while True:
            if j == len(word) + 1 or index_in_alphabet(word[j - 1], list_alphabet) < index_in_alphabet(word[i - 1],
                                                                                                       list_alphabet):
                while k < i:
                    # print(word[k:k + j - i])
                    w = word[k:k + j - i]
                    if len(w) <= C:
                        CFL_list.append(word[k:k + j - i])
                    else:
                        ICFL_list_recursive = ICFL_recursive(w)

                        # Insert << to indicate the begin of the subdecomposition of w
                        CFL_list.append("<<")
                        for v in ICFL_list_recursive:
                            CFL_list.append(v)
                        # Insert >> to indicate the end of the subdecomposition of w
                        CFL_list.append(">>")

                    k = k + j - i
                break
            else:
                # if word[j-1] > word[i-1]:
                if index_in_alphabet(word[j - 1], list_alphabet) > index_in_alphabet(word[i - 1], list_alphabet):
                    i = k + 1
                else:
                    i = i + 1
                j = j + 1

    return CFL_list


# ---------------------------------------------------------------------------------------------



# ------------------------ CFL_icfl_cfl ---------------------------------------------------------------------
# CFL factorization - ICFL subdecomposition - CFL factorization
def CFL_icfl_cfl(word, C):
    """
    CFL Duval's algorithm.
    """
    CFL_list = []
    k = 0
    while k < len(word):
        i = k + 1
        j = k + 2
        while True:
            if j == len(word) + 1 or word[j - 1] < word[i - 1]:
                while k < i:
                    # print(word[k:k + j - i])
                    w = word[k:k + j - i]
                    if len(w) <= C:
                        CFL_list.append(word[k:k + j - i])
                    else:
                        ICFL_list_recursive = ICFL_cfl(w, C)

                        # Insert << to indicate the begin of the subdecomposition of w
                        CFL_list.append("<<")
                        for v in ICFL_list_recursive:
                            CFL_list.append(v)
                        # Insert >> to indicate the end of the subdecomposition of w
                        CFL_list.append(">>")

                    k = k + j - i
                break
            else:
                # if word[j-1] > word[i-1]:
                if word[j - 1] > word[i - 1]:
                    i = k + 1
                else:
                    i = i + 1
                j = j + 1

    return CFL_list


# CFL factorization - ICFL subdecomposition - on a specific algorithm
def CFL_icfl_cfl_for_alphabet(word, C, list_alphabet):
    """
    CFL Duval's algorithm.
    """
    CFL_list = []
    k = 0
    while k < len(word):
        i = k + 1
        j = k + 2
        while True:
            if j == len(word) + 1 or index_in_alphabet(word[j - 1], list_alphabet) < index_in_alphabet(word[i - 1],
                                                                                                       list_alphabet):
                while k < i:
                    # print(word[k:k + j - i])
                    w = word[k:k + j - i]
                    if len(w) <= C:
                        CFL_list.append(word[k:k + j - i])
                    else:
                        ICFL_list_recursive = ICFL_cfl_for_alphabet(w, C, list_alphabet)

                        # Insert << to indicate the begin of the subdecomposition of w
                        CFL_list.append("<<")
                        for v in ICFL_list_recursive:
                            CFL_list.append(v)
                        # Insert >> to indicate the end of the subdecomposition of w
                        CFL_list.append(">>")

                    k = k + j - i
                break
            else:
                # if word[j-1] > word[i-1]:
                if index_in_alphabet(word[j - 1], list_alphabet) > index_in_alphabet(word[i - 1], list_alphabet):
                    i = k + 1
                else:
                    i = i + 1
                j = j + 1

    return CFL_list


# ---------------------------------------------------------------------------------------------

# ICFL recursive factorization - CFL subdecomposition
def ICFL_cfl(word, C):
    """In this version of ICFL, we don't execute compute_br - one only O(n) scanning of word"""
    #br_list = []
    icfl_list = []

    compute_icfl_recursive(word, icfl_list)

    ICFL_cfl_list = []
    for w in icfl_list:
        if len(w) <= C:
            ICFL_cfl_list.append(w)
        else:
            CFL_list = CFL(w)

            # Insert << to indicate the begin of the subdecomposition of w
            ICFL_cfl_list.append("<<")
            for v in CFL_list:
                ICFL_cfl_list.append(v)
            # Insert >> to indicate the end of the subdecomposition of w
            ICFL_cfl_list.append(">>")

    return ICFL_cfl_list


# ICFL recursive factorization - CFL subdecomposition - for a specific alphabet
def ICFL_cfl_for_alphabet(word, C, list_alphabet):
    """In this version of ICFL, we don't execute compute_br - one only O(n) scanning of word"""
    br_list = []
    icfl_list = []

    compute_icfl_recursive_for_alphabet(word, icfl_list, list_alphabet)

    ICFL_cfl_list = []
    for w in icfl_list:
        if len(w) <= C:
            ICFL_cfl_list.append(w)
        else:
            CFL_list = CFL_for_alphabet(w, list_alphabet)

            # Insert << to indicate the begin of the subdecomposition of w
            ICFL_cfl_list.append("<<")
            for v in CFL_list:
                ICFL_cfl_list.append(v)
            # Insert >> to indicate the end of the subdecomposition of w
            ICFL_cfl_list.append(">>")

    return ICFL_cfl_list


# ----------------------------------------------------------------------------------------------------------------------------------------------------


# ICFL recursive factorization - CFL subdecomposition - ICFL recursive factorization
def ICFL_cfl_icfl(word, C):
    """In this version of ICFL, we don't execute compute_br - one only O(n) scanning of word"""
    br_list = []
    icfl_list = []

    compute_icfl_cfl_icfl(word, br_list, icfl_list)

    ICFL_cfl_list = []
    for w in icfl_list:
        if len(w) <= C:
            ICFL_cfl_list.append(w)
        else:
            CFL_list = CFL_icfl(w, C)

            # Insert << to indicate the begin of the subdecomposition of w
            ICFL_cfl_list.append("<<")
            for v in CFL_list:
                ICFL_cfl_list.append(v)
            # Insert >> to indicate the end of the subdecomposition of w
            ICFL_cfl_list.append(">>")

    return ICFL_cfl_list


def compute_icfl_cfl_icfl(word, br_list, icfl_list):
    # At each step compute the current bre
    pre_pair = find_pre(word)
    current_bre_quad = find_bre(pre_pair)
    br_list.append(current_bre_quad)

    if current_bre_quad[1] == '' and current_bre_quad[0].find('$') >= 0:
        w = current_bre_quad[0]
        icfl_list.insert(0, w[:len(w) - 1])
        return
    else:
        compute_icfl_recursive(current_bre_quad[1] + current_bre_quad[2], br_list, icfl_list)
        if len(icfl_list[0]) > current_bre_quad[3]:
            icfl_list.insert(0, current_bre_quad[0])
        else:
            icfl_list[0] = current_bre_quad[0] + icfl_list[0]
        return


# ICFL recursive factorization - CFL subdecomposition - ICFL recursive factorization - for a specific alphabet
def ICFL_cfl_icfl_for_alphabet(word, C, list_alphabet):
    """In this version of ICFL, we don't execute compute_br - one only O(n) scanning of word"""
    br_list = []
    icfl_list = []

    compute_icfl_cfl_icfl_for_alphabet(word, br_list, icfl_list, list_alphabet)

    ICFL_cfl_list = []
    for w in icfl_list:
        if len(w) <= C:
            ICFL_cfl_list.append(w)
        else:
            CFL_list = CFL_icfl_for_alphabet(w, C, list_alphabet)

            # Insert << to indicate the begin of the subdecomposition of w
            ICFL_cfl_list.append("<<")
            for v in CFL_list:
                ICFL_cfl_list.append(v)
            # Insert >> to indicate the end of the subdecomposition of w
            ICFL_cfl_list.append(">>")

    return ICFL_cfl_list


def compute_icfl_cfl_icfl_for_alphabet(word, br_list, icfl_list, list_alphabet):
    # At each step compute the current bre
    pre_pair = find_pre_for_alphabet(word, list_alphabet)
    current_bre_quad = find_bre_for_alphabet(pre_pair, list_alphabet)
    br_list.append(current_bre_quad)

    if current_bre_quad[1] == '' and current_bre_quad[0].find('$') >= 0:
        w = current_bre_quad[0]
        icfl_list.insert(0, w[:len(w) - 1])
        return
    else:
        compute_icfl_recursive(current_bre_quad[1] + current_bre_quad[2], br_list, icfl_list)
        if len(icfl_list[0]) > current_bre_quad[3]:
            icfl_list.insert(0, current_bre_quad[0])
        else:
            icfl_list[0] = current_bre_quad[0] + icfl_list[0]
        return


# --------------------------------------------------------------------------------------------------------------



# ------------------------ CFL_cflin ---------------------------------------------------------------------
# CFL factorization - CFL_in subdecomposition
def CFL_cflin(word, C):
    """
    CFL Duval's algorithm.
    """
    CFL_list = []
    k = 0
    while k < len(word):
        i = k + 1
        j = k + 2
        while True:
            if j == len(word) + 1 or word[j - 1] < word[i - 1]:
                while k < i:
                    # print(word[k:k + j - i])
                    w = word[k:k + j - i]
                    if len(w) <= C:
                        CFL_list.append(word[k:k + j - i])
                    else:
                        CFL_in_list = CFL_for_alphabet(word, ['T', 'N', 'G', 'C', 'A'])

                        # Insert << to indicate the begin of the subdecomposition of w
                        CFL_list.append("<<")
                        for v in CFL_in_list:
                            CFL_list.append(v)
                        # Insert >> to indicate the end of the subdecomposition of w
                        CFL_list.append(">>")

                    k = k + j - i
                break
            else:
                # if word[j-1] > word[i-1]:
                if word[j - 1] > word[i - 1]:
                    i = k + 1
                else:
                    i = i + 1
                j = j + 1

    return CFL_list


# CFL factorization - CFL_in subdecomposition - on a specific algorithm
def CFL_cflin_for_alphabet(word, C, list_alphabet):
    """
    CFL Duval's algorithm.
    """
    CFL_list = []
    k = 0
    while k < len(word):
        i = k + 1
        j = k + 2
        while True:
            if j == len(word) + 1 or index_in_alphabet(word[j - 1], list_alphabet) < index_in_alphabet(word[i - 1],
                                                                                                       list_alphabet):
                while k < i:
                    # print(word[k:k + j - i])
                    w = word[k:k + j - i]
                    if len(w) <= C:
                        CFL_list.append(word[k:k + j - i])
                    else:
                        CFL_in_list = CFL_for_alphabet(word, list_alphabet[::-1])

                        # Insert << to indicate the begin of the subdecomposition of w
                        CFL_list.append("<<")
                        for v in CFL_in_list:
                            CFL_list.append(v)
                        # Insert >> to indicate the end of the subdecomposition of w
                        CFL_list.append(">>")

                    k = k + j - i
                break
            else:
                # if word[j- 1] > word[i-1]:
                if index_in_alphabet(word[j - 1], list_alphabet) > index_in_alphabet(word[i - 1], list_alphabet):
                    i = k + 1
                else:
                    i = i + 1
                j = j + 1

    return CFL_list


# ---------------------------------------------------------------------------------------------



# ------------------------ CFL_cflin_cfl ---------------------------------------------------------------------
# CFL factorization - CFL_in subdecomposition - CFL factorization
def CFL_cflin_cfl(word, C):
    """
    CFL Duval's algorithm.
    """
    CFL_list = []
    k = 0
    while k < len(word):
        i = k + 1
        j = k + 2
        while True:
            if j == len(word) + 1 or word[j - 1] < word[i - 1]:
                while k < i:
                    # print(word[k:k + j - i])
                    w = word[k:k + j - i]
                    if len(w) <= C:
                        CFL_list.append(word[k:k + j - i])
                    else:
                        CFL_in_list = CFL_cflin(w, C)

                        # Insert << to indicate the begin of the subdecomposition of w
                        CFL_list.append("<<")
                        for v in CFL_in_list:
                            CFL_list.append(v)
                        # Insert >> to indicate the end of the subdecomposition of w
                        CFL_list.append(">>")

                    k = k + j - i
                break
            else:
                # if word[j-1] > word[i-1]:
                if word[j - 1] > word[i - 1]:
                    i = k + 1
                else:
                    i = i + 1
                j = j + 1

    return CFL_list


# CFL factorization - CFL_in subdecomposition - on a specific algorithm
def CFL_cflin_cfl_for_alphabet(word, C, list_alphabet):
    """
    CFL Duval's algorithm.
    """
    CFL_list = []
    k = 0
    while k < len(word):
        i = k + 1
        j = k + 2
        while True:
            if j == len(word) + 1 or index_in_alphabet(word[j - 1], list_alphabet) < index_in_alphabet(word[i - 1],
                                                                                                       list_alphabet):
                while k < i:
                    # print(word[k:k + j - i])
                    w = word[k:k + j - i]
                    if len(w) <= C:
                        CFL_list.append(word[k:k + j - i])
                    else:
                        CFL_in_list = CFL_cflin_for_alphabet(w, C, list_alphabet)

                        # Insert << to indicate the begin of the subdecomposition of w
                        CFL_list.append("<<")
                        for v in CFL_in_list:
                            CFL_list.append(v)
                        # Insert >> to indicate the end of the subdecomposition of w
                        CFL_list.append(">>")

                    k = k + j - i
                break
            else:
                # if word[j-1] > word[i-1]:
                if index_in_alphabet(word[j - 1], list_alphabet) > index_in_alphabet(word[i - 1], list_alphabet):
                    i = k + 1
                else:
                    i = i + 1
                j = j + 1

    return CFL_list
