import time

a = 1
b = 2
c = a + b


# print(c)


def testFunction():
    for i in range(1000):
        for g in range(5):
            time.sleep(1)
            print(i)


print('HI')

testFunction()
