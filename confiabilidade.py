def serie(R):
    """ Returns the Confiability of components in series """
    ans= 1
    for i in R:
        ans=ans*i

    return (ans)

def paralelo(R):
    """ Returns the Total Confiability of components in parallel"""
    mult = 1
    for i in R:
        mult = mult * (1-i)
    8
    return (1-mult)

