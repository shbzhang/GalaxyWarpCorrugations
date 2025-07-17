import numpy as np

### P1
i = 9

while True:
	if i % 9 == 0:
		ri = int(str(i)[::-1])
		if ri % 7 == 0:
			break
	i+=1

print(i, ri)


### P2
def weird_append(val, lst=[]):
    lst.append(val)
    return lst

print(weird_append(1))   # Output?
print(weird_append(2))   # Output?
print(weird_append(3))   # Output?




### P3
import random

def unfair_shuffle(lst):
    return sorted(lst, key=lambda x: random.random())

s = np.zeros(10)
for i in range(1000):
	lst = list(range(10))
	shuffled = unfair_shuffle(lst)
	s += shuffled

print(s)



### P4
def find_missing_number(a):
	s = sum(a)
	return sum(range(1,101))-s


import random
def random_missing():
	lst = list(range(1,101))
	drop = random.randint(1,100)
	lst.remove(drop)
	random.shuffle(lst)
	return lst

a = random_missing()
print(find_missing_number(a))
print(set(range(1,101)) - set(a))



### P5
def strange(n):
    if n == 0:
        return 0
    return n + strange(n - 5)

print(strange(5))  # What's the output and why?



### P6
mylist = ['apple', 'banana', 'cherry']
x = frozenset(mylist)
x[1] = "strawberry"


a = {1: 'one', 2: 'two'}
b = frozenset(a.items())
c = dict(b)
print(a,b,c)
