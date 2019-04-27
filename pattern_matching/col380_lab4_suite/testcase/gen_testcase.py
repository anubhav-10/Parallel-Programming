#!/usr/bin/python3

#########################################################################
# Generate the input file containing text and pattern in required       #
# format and save to file named testcase_<n>_<num_patterns>             #
#                                                                       #
# Parameters:                                                           #
#   n               :length of text                                     #
#   num_patterns    :number of patterns to be searched                  #
#   min_p           :minimum period length                              #
#   min_m           :minimum period length                              #
# Format of output file:                                                #
#   -----------------------------------------------------------------   #
#   | n num_patterns                                                    #
#   | text                                                              #
#   | m[0] m[1] m[2] ... m[num_patterns-1]                              #
#   | p[0] p[1] p[2] ... p[num_patterns-1]                              #
#   | pattern[0]                                                        #
#   | pattern[1]                                                        #
#   | ...                                                               #
#   | pattern[num_pattern-1]                                            #
#   -----------------------------------------------------------------   #
#                                                                       #
#########################################################################

import random
import string
import re

# adjust these parameters to generate testcases
n            = 10000    # length of text
num_patterns = 10       # number of patterns to be searched
min_p        = 2        # minimum period length
min_m        = 5        # minimum pattern length
num_chars	 = 26		# number of different characters in text and patterns,
						# can't be more than 26
max_m		 = 100		# maximum length of patterns, can be upto n

def randomString(stringLength):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase[:num_chars]
    return ''.join(random.choice(letters) for i in range(stringLength))

def isNonPeriodic(testString):
	x = re.search(r'\b(\w*)(\w+\1)\2+\b', testString)
	if x:
		print(testString, 'is periodic')
		return True
	else:
		print(testString, 'is not periodic')
		return False


m_set = []
p_set = []
pattern_set = []

# generate length of patterns
for i in range(num_patterns):
	m = random.randint(min_m, max_m)
	m_set.append(m)

# generate period length of patterns
for i in range(num_patterns):
	p = random.randint(min_p, m_set[i]//2)
	p_set.append(p)

# generate patterns
for i in range(num_patterns):
	period_str = randomString(p_set[i])
	while isNonPeriodic(period_str):
		period_str = randomString(p_set[i])
	pattern = period_str * (m_set[i] // p_set[i])
	pattern = pattern + period_str[0 : (m_set[i]-len(pattern))]
	pattern_set.append(pattern)

# generate text
text = randomString(n)
max_matches = n//max_m
for i in range (max_matches):
	pattern_no = random.randint(0, num_patterns-1)
	pos = random.randint(0, n - m_set[pattern_no])
	text = text[0:pos] + pattern_set[pattern_no] + text[pos+m_set[pattern_no]: ]

filename = 'testcase_' + str(n) + '_' + str(num_patterns)     #output filename
file = open(filename, 'w')

# write size of length of text and number of patterns in first line of file
file.write(str(n) + ' ' +str(num_patterns) + '\n')

# write text
file.write(text + '\n')

# write length of patterns
for i in range(num_patterns):
	file.write(str(m_set[i]) + ' ')
file.write('\n')

# generate and write period length of patterns
for i in range(num_patterns):
	file.write(str(p_set[i]) + ' ')
file.write('\n')

# write patterns
for i in range(num_patterns):
	file.write(pattern_set[i] + '\n')

file.close()
