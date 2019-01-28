def nextregular(n):
    while not checksize(n): n+=1
    return n

def checksize(n):
    while not (n%16): n/=16
    while not (n%13): n/=13
    while not (n%11): n/=11
    while not (n%9): n/=9
    while not (n%7): n/=7
    while not (n%5): n/=5
    while not (n%3): n/=3
    while not (n%2): n/=2
    return (1 if n == 1 else 0)
