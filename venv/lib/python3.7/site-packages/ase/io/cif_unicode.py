'''
Conversion of text from a Crystallographic Information File (CIF) format to
unicode. CIF text is neither unicode nor bibtex/latex code.

Rules for character formatting in CIF files are specified at:
https://www.iucr.org/resources/cif/spec/version1.1/semantics
'''

import re
import html


subs_dict = {
    '\r': '',            # Windows line ending
    '\t': ' ',           # tabs

    r'\a': u'\u03b1',    # alpha
    r'\b': u'\u03b2',    # beta
    r'\g': u'\u03b3',    # gamma
    r'\d': u'\u03b4',    # delta
    r'\e': u'\u03b5',    # epsilon
    r'\z': u'\u03b6',    # zeta
    r'\h': u'\u03b7',    # eta
    r'\q': u'\u03b8',    # theta
    r'\i': u'\u03b9',    # iota
    r'\k': u'\u03ba',    # kappa
    r'\l': u'\u03bb',    # lambda
    r'\m': u'\u03bc',    # mu
    r'\n': u'\u03bd',    # nu
    r'\x': u'\u03be',    # xi
    r'\o': u'\u03bf',    # omicron
    r'\p': u'\u03c0',    # pi
    r'\r': u'\u03c1',    # rho
    r'\s': u'\u03c3',    # sigma
    r'\t': u'\u03c4',    # tau
    r'\u': u'\u03c5',    # upsilon
    r'\f': u'\u03c6',    # phi
    r'\c': u'\u03c7',    # chi
    r'\y': u'\u03c8',    # psi
    r'\w': u'\u03c9',    # omega
    r'\A': u'\u0391',    # Alpha
    r'\B': u'\u0392',    # Beta
    r'\G': u'\u0393',    # Gamma
    r'\D': u'\u0394',    # Delta
    r'\E': u'\u0395',    # Epsilon
    r'\Z': u'\u0396',    # Zeta
    r'\H': u'\u0397',    # Eta
    r'\Q': u'\u0398',    # Theta
    r'\I': u'\u0399',    # Ioto
    r'\K': u'\u039a',    # Kappa
    r'\L': u'\u039b',    # Lambda
    r'\M': u'\u039c',    # Mu
    r'\N': u'\u039d',    # Nu
    r'\X': u'\u039e',    # Xi
    r'\O': u'\u039f',    # Omicron
    r'\P': u'\u03a0',    # Pi
    r'\R': u'\u03a1',    # Rho
    r'\S': u'\u03a3',    # Sigma
    r'\T': u'\u03a4',    # Tau
    r'\U': u'\u03a5',    # Upsilon
    r'\F': u'\u03a6',    # Phi
    r'\C': u'\u03a7',    # Chi
    r'\Y': u'\u03a8',    # Psi
    r'\W': u'\u03a9',    # Omega

    r'\%a': u'\u00e5',   # a-ring
    r'\/o': u'\u00f8',   # o-slash
    r'\?i': u'\u0131',   # dotless i
    r'\/l': u'\u0142',   # Polish l
    r'\&s': u'\u00df',   # German eszett
    r'\/d': u'\u0111',   # barred d

    r'\%A': u'\u00c5',   # A-ring
    r'\/O': u'\u00d8',   # O-slash
    r'\?I': 'I',         # dotless I
    r'\/L': u'\u0141',   # Polish L
    r'\&S': u'\u1e9e',   # German Eszett
    r'\/D': u'\u0110',   # barred D

    r'\%': u'\u00b0',           # degree
    r'--': u'\u2013',           # dash
    r'---': u'\u2014',          # single bond
    r'\\db': u'\u003d',         # double bond
    r'\\tb': u'\u2261',         # triple bond
    r'\\ddb': u'\u2248',        # delocalized double bond
    r'\\sim': '~',
    r'\\simeq': u'\u2243',
    r'\\infty': u'\u221e',      # infinity

    r'\\times': u'\u00d7',
    r'+-': u'\u00b1',           # plusminus
    r'-+': u'\u2213',           # minusplus
    r'\\square': u'\u25a0',
    r'\\neq': u'\u2660',
    r'\\rangle': u'\u3009',
    r'\\langle': u'\u3008',
    r'\\rightarrow': u'\u2192',
    r'\\leftarrow': u'\u2190',

    r"\'A": u'\u00c1',  # A acute
    r"\'E": u'\u00c9',  # E acute
    r"\'I": u'\u00cd',  # I acute
    r"\'O": u'\u00d3',  # O acute
    r"\'U": u'\u00da',  # U acute
    r"\'Y": u'\u00dd',  # Y acute
    r"\'a": u'\u00e1',  # a acute
    r"\'e": u'\u00e9',  # e acute
    r"\'i": u'\u00ed',  # i acute
    r"\'o": u'\u00f3',  # o acute
    r"\'u": u'\u00fa',  # u acute
    r"\'y": u'\u00fd',  # y acute
    r"\'C": u'\u0106',  # C acute
    r"\'c": u'\u0107',  # c acute
    r"\'L": u'\u0139',  # L acute
    r"\'l": u'\u013a',  # l acute
    r"\'N": u'\u0143',  # N acute
    r"\'n": u'\u0144',  # n acute
    r"\'R": u'\u0154',  # R acute
    r"\'r": u'\u0155',  # r acute
    r"\'S": u'\u015a',  # S acute
    r"\'s": u'\u015b',  # s acute
    r"\'Z": u'\u0179',  # Z acute
    r"\'z": u'\u017a',  # z acute
    r"\'G": u'\u01f4',  # G acute
    r"\'g": u'\u01f5',  # g acute
    r"\'K": u'\u1e30',  # K acute
    r"\'k": u'\u1e31',  # k acute
    r"\'M": u'\u1e3e',  # M acute
    r"\'m": u'\u1e3f',  # m acute
    r"\'P": u'\u1e54',  # P acute
    r"\'p": u'\u1e55',  # p acute
    r"\'W": u'\u1e82',  # W acute
    r"\'w": u'\u1e83',  # w acute
    r'\;A': u'\u0104',  # A ogonek
    r'\;a': u'\u0105',  # a ogonek
    r'\;E': u'\u0118',  # E ogonek
    r'\;e': u'\u0119',  # e ogonek
    r'\;I': u'\u012e',  # I ogonek
    r'\;i': u'\u012f',  # i ogonek
    r'\;U': u'\u0172',  # U ogonek
    r'\;u': u'\u0173',  # u ogonek
    r'\;O': u'\u01ea',  # O ogonek
    r'\;o': u'\u01eb',  # o ogonek
    r'\.C': u'\u010a',  # C dot above
    r'\.c': u'\u010b',  # c dot above
    r'\.E': u'\u0116',  # E dot above
    r'\.e': u'\u0117',  # e dot above
    r'\.G': u'\u0120',  # G dot above
    r'\.g': u'\u0121',  # g dot above
    r'\.I': u'\u0130',  # I dot above
    r'\.Z': u'\u017b',  # Z dot above
    r'\.z': u'\u017c',  # z dot above
    r'\.A': u'\u0226',  # A dot above
    r'\.a': u'\u0227',  # a dot above
    r'\.O': u'\u022e',  # O dot above
    r'\.o': u'\u022f',  # o dot above
    r'\.B': u'\u1e02',  # B dot above
    r'\.b': u'\u1e03',  # b dot above
    r'\.D': u'\u1e0a',  # D dot above
    r'\.d': u'\u1e0b',  # d dot above
    r'\.F': u'\u1e1e',  # F dot above
    r'\.f': u'\u1e1f',  # f dot above
    r'\.H': u'\u1e22',  # H dot above
    r'\.h': u'\u1e23',  # h dot above
    r'\.M': u'\u1e40',  # M dot above
    r'\.m': u'\u1e41',  # m dot above
    r'\.N': u'\u1e44',  # N dot above
    r'\.n': u'\u1e45',  # n dot above
    r'\.P': u'\u1e56',  # P dot above
    r'\.p': u'\u1e57',  # p dot above
    r'\.R': u'\u1e58',  # R dot above
    r'\.r': u'\u1e59',  # r dot above
    r'\.S': u'\u1e60',  # S dot above
    r'\.s': u'\u1e61',  # s dot above
    r'\.T': u'\u1e6a',  # T dot above
    r'\.t': u'\u1e6b',  # t dot above
    r'\.W': u'\u1e86',  # W dot above
    r'\.w': u'\u1e87',  # w dot above
    r'\.X': u'\u1e8a',  # X dot above
    r'\.x': u'\u1e8b',  # x dot above
    r'\.Y': u'\u1e8e',  # Y dot above
    r'\.y': u'\u1e8f',  # y dot above
    r'\(A': u'\u0102',  # A breve
    r'\(a': u'\u0103',  # a breve
    r'\(E': u'\u0114',  # E breve
    r'\(e': u'\u0115',  # e breve
    r'\(G': u'\u011e',  # G breve
    r'\(g': u'\u011f',  # g breve
    r'\(I': u'\u012c',  # I breve
    r'\(i': u'\u012d',  # i breve
    r'\(O': u'\u014e',  # O breve
    r'\(o': u'\u014f',  # o breve
    r'\(U': u'\u016c',  # U breve
    r'\(u': u'\u016d',  # u breve
    r'\=A': u'\u0100',  # A macron
    r'\=a': u'\u0101',  # a macron
    r'\=E': u'\u0112',  # E macron
    r'\=e': u'\u0113',  # e macron
    r'\=I': u'\u012a',  # I macron
    r'\=i': u'\u012b',  # i macron
    r'\=O': u'\u014c',  # O macron
    r'\=o': u'\u014d',  # o macron
    r'\=U': u'\u016a',  # U macron
    r'\=u': u'\u016b',  # u macron
    r'\=Y': u'\u0232',  # Y macron
    r'\=y': u'\u0233',  # y macron
    r'\=G': u'\u1e20',  # G macron
    r'\=g': u'\u1e21',  # g macron
    r'\^A': u'\u00c2',  # A circumflex
    r'\^E': u'\u00ca',  # E circumflex
    r'\^I': u'\u00ce',  # I circumflex
    r'\^O': u'\u00d4',  # O circumflex
    r'\^U': u'\u00db',  # U circumflex
    r'\^a': u'\u00e2',  # a circumflex
    r'\^e': u'\u00ea',  # e circumflex
    r'\^i': u'\u00ee',  # i circumflex
    r'\^o': u'\u00f4',  # o circumflex
    r'\^u': u'\u00fb',  # u circumflex
    r'\^C': u'\u0108',  # C circumflex
    r'\^c': u'\u0109',  # c circumflex
    r'\^G': u'\u011c',  # G circumflex
    r'\^g': u'\u011d',  # g circumflex
    r'\^H': u'\u0124',  # H circumflex
    r'\^h': u'\u0125',  # h circumflex
    r'\^J': u'\u0134',  # J circumflex
    r'\^j': u'\u0135',  # j circumflex
    r'\^S': u'\u015c',  # S circumflex
    r'\^s': u'\u015d',  # s circumflex
    r'\^W': u'\u0174',  # W circumflex
    r'\^w': u'\u0175',  # w circumflex
    r'\^Y': u'\u0176',  # Y circumflex
    r'\^y': u'\u0177',  # y circumflex
    r'\^Z': u'\u1e90',  # Z circumflex
    r'\^z': u'\u1e91',  # z circumflex
    r'\"A': u'\u00c4',  # A diaeresis
    r'\"E': u'\u00cb',  # E diaeresis
    r'\"I': u'\u00cf',  # I diaeresis
    r'\"O': u'\u00d6',  # O diaeresis
    r'\"U': u'\u00dc',  # U diaeresis
    r'\"a': u'\u00e4',  # a diaeresis
    r'\"e': u'\u00eb',  # e diaeresis
    r'\"i': u'\u00ef',  # i diaeresis
    r'\"o': u'\u00f6',  # o diaeresis
    r'\"u': u'\u00fc',  # u diaeresis
    r'\"y': u'\u00ff',  # y diaeresis
    r'\"Y': u'\u0178',  # Y diaeresis
    r'\"H': u'\u1e26',  # H diaeresis
    r'\"h': u'\u1e27',  # h diaeresis
    r'\"W': u'\u1e84',  # W diaeresis
    r'\"w': u'\u1e85',  # w diaeresis
    r'\"X': u'\u1e8c',  # X diaeresis
    r'\"x': u'\u1e8d',  # x diaeresis
    r'\"t': u'\u1e97',  # t diaeresis
    r'\~A': u'\u00c3',  # A tilde
    r'\~N': u'\u00d1',  # N tilde
    r'\~O': u'\u00d5',  # O tilde
    r'\~a': u'\u00e3',  # a tilde
    r'\~n': u'\u00f1',  # n tilde
    r'\~o': u'\u00f5',  # o tilde
    r'\~I': u'\u0128',  # I tilde
    r'\~i': u'\u0129',  # i tilde
    r'\~U': u'\u0168',  # U tilde
    r'\~u': u'\u0169',  # u tilde
    r'\~V': u'\u1e7c',  # V tilde
    r'\~v': u'\u1e7d',  # v tilde
    r'\~E': u'\u1ebc',  # E tilde
    r'\~e': u'\u1ebd',  # e tilde
    r'\~Y': u'\u1ef8',  # Y tilde
    r'\~y': u'\u1ef9',  # y tilde
    r'\<C': u'\u010c',  # C caron
    r'\<c': u'\u010d',  # c caron
    r'\<D': u'\u010e',  # D caron
    r'\<d': u'\u010f',  # d caron
    r'\<E': u'\u011a',  # E caron
    r'\<e': u'\u011b',  # e caron
    r'\<L': u'\u013d',  # L caron
    r'\<l': u'\u013e',  # l caron
    r'\<N': u'\u0147',  # N caron
    r'\<n': u'\u0148',  # n caron
    r'\<R': u'\u0158',  # R caron
    r'\<r': u'\u0159',  # r caron
    r'\<S': u'\u0160',  # S caron
    r'\<s': u'\u0161',  # s caron
    r'\<T': u'\u0164',  # T caron
    r'\<t': u'\u0165',  # t caron
    r'\<Z': u'\u017d',  # Z caron
    r'\<z': u'\u017e',  # z caron
    r'\<A': u'\u01cd',  # A caron
    r'\<a': u'\u01ce',  # a caron
    r'\<I': u'\u01cf',  # I caron
    r'\<i': u'\u01d0',  # i caron
    r'\<O': u'\u01d1',  # O caron
    r'\<o': u'\u01d2',  # o caron
    r'\<U': u'\u01d3',  # U caron
    r'\<u': u'\u01d4',  # u caron
    r'\<G': u'\u01e6',  # G caron
    r'\<g': u'\u01e7',  # g caron
    r'\<K': u'\u01e8',  # K caron
    r'\<k': u'\u01e9',  # k caron
    r'\<j': u'\u01f0',  # j caron
    r'\<H': u'\u021e',  # H caron
    r'\<h': u'\u021f',  # h caron
    r'\>O': u'\u0150',  # O double acute
    r'\>o': u'\u0151',  # o double acute
    r'\>U': u'\u0170',  # U double acute
    r'\>u': u'\u0171',  # u double acute
    r'\,C': u'\u00c7',  # C cedilla
    r'\,c': u'\u00e7',  # c cedilla
    r'\,G': u'\u0122',  # G cedilla
    r'\,g': u'\u0123',  # g cedilla
    r'\,K': u'\u0136',  # K cedilla
    r'\,k': u'\u0137',  # k cedilla
    r'\,L': u'\u013b',  # L cedilla
    r'\,l': u'\u013c',  # l cedilla
    r'\,N': u'\u0145',  # N cedilla
    r'\,n': u'\u0146',  # n cedilla
    r'\,R': u'\u0156',  # R cedilla
    r'\,r': u'\u0157',  # r cedilla
    r'\,S': u'\u015e',  # S cedilla
    r'\,s': u'\u015f',  # s cedilla
    r'\,T': u'\u0162',  # T cedilla
    r'\,t': u'\u0163',  # t cedilla
    r'\,E': u'\u0228',  # E cedilla
    r'\,e': u'\u0229',  # e cedilla
    r'\,D': u'\u1e10',  # D cedilla
    r'\,d': u'\u1e11',  # d cedilla
    r'\,H': u'\u1e28',  # H cedilla
    r'\,h': u'\u1e29',  # h cedilla
    r'\`A': u'\u00c0',  # A grave
    r'\`E': u'\u00c8',  # E grave
    r'\`I': u'\u00cc',  # I grave
    r'\`O': u'\u00d2',  # O grave
    r'\`U': u'\u00d9',  # U grave
    r'\`a': u'\u00e0',  # a grave
    r'\`e': u'\u00e8',  # e grave
    r'\`i': u'\u00ec',  # i grave
    r'\`o': u'\u00f2',  # o grave
    r'\`u': u'\u00f9',  # u grave
    r'\`N': u'\u01f8',  # N grave
    r'\`n': u'\u01f9',  # n grave
    r'\`W': u'\u1e80',  # W grave
    r'\`w': u'\u1e81',  # w grave
    r'\`Y': u'\u1ef2',  # Y grave
    r'\`y': u'\u1ef3',  # y grave
}

superscript_dict = {
    '0': u'\u2070',  # superscript 0
    '1': u'\u00b9',  # superscript 1
    '2': u'\u00b2',  # superscript 2
    '3': u'\u00b3',  # superscript 3
    '4': u'\u2074',  # superscript 4
    '5': u'\u2075',  # superscript 5
    '6': u'\u2076',  # superscript 6
    '7': u'\u2077',  # superscript 7
    '8': u'\u2078',  # superscript 8
    '9': u'\u2079',  # superscript 9
}

subscript_dict = {
    '0': u'\u2080',  # subscript 0
    '1': u'\u2081',  # subscript 1
    '2': u'\u2082',  # subscript 2
    '3': u'\u2083',  # subscript 3
    '4': u'\u2084',  # subscript 4
    '5': u'\u2085',  # subscript 5
    '6': u'\u2086',  # subscript 6
    '7': u'\u2087',  # subscript 7
    '8': u'\u2088',  # subscript 8
    '9': u'\u2089',  # subscript 9
}


def replace_subscript(s, subscript=True):

    target = '~'
    rdict = subscript_dict
    if not subscript:
        target = '^'
        rdict = superscript_dict

    replaced = []
    inside = False
    for char in s:
        if char == target:
            inside = not inside
        elif not inside:
            replaced += [char]
        # note: do not use char.isdigit - this also matches (sub/super)scripts
        elif char in rdict:
            replaced += [rdict[char]]
        else:
            replaced += [char]

    return ''.join(replaced)


def multiple_replace(text, adict):
    rx = re.compile('|'.join(map(re.escape, adict)))

    def one_xlat(match):
        return adict[match.group(0)]

    return rx.sub(one_xlat, text)


def format_unicode(s):
    """Converts a string in CIF text-format to unicode.  Any HTML tags
    contained in the string are removed.  HTML numeric character references
    are unescaped (i.e. converted to unicode).

    Parameters:

    s: string
        The CIF text string to convert

    Returns:

    u: string
        A unicode formatted string.
    """

    s = html.unescape(s)
    s = multiple_replace(s, subs_dict)
    tagclean = re.compile('<.*?>')
    return re.sub(tagclean, '', s)


def handle_subscripts(s):
    s = replace_subscript(s, subscript=True)
    s = replace_subscript(s, subscript=False)
    return s
