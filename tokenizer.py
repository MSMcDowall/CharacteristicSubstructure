import graphADT
import re

#probelem getting the regular expression to check the aromatic set first and then [A-Z][a-z]?
# |\d|b|c|n|o|p|s|se|as|[A-Z][a-z]?
delimiter_pattern = re.compile('([-=#\$:/\\ \(\).\[\]|\w])')
broken = delimiter_pattern.findall('COcc(c.c#c)2cCNTi')

# concatenate capital letter with small letter(not aromatic) to make atomic symbol
aromatic = re.compile('b|c|n|o|p|s|se|as')
lowercase = re.compile('[a-z]')

new_broken = []
for token in broken:
    if re.match(lowercase, token) and not re.match(aromatic, token):
        position_lc = broken.index(token)
        new_broken.append(broken[position_lc-1]+broken[position_lc])
    else:
        new_broken.append(token)


print broken
print new_broken
