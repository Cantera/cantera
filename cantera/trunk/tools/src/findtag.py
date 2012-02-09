import sys

def cvstag(dir):
    path = dir+"/CVS/Tag"
    try:
        f = open(path,'r')
        if f:
            s = f.readlines()
            if s:
                return s[0][1:-1]
            else:
                return ""
        else:
            return ""
    except:
        return ""

if __name__ == "__main__":
    dir = sys.argv[1]
    print cvstag(dir)
    
