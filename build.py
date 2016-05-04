#!/usr/bin/env python

import sys, os



def grab(srcdir, line):
    #print line,
    line = line.split("{")[1]
    name = line.split("}")[0]
    print ("cp %s/%s %s"%(srcdir, name, name))
    os.system("cp %s/%s %s"%(srcdir, name, name))


def process(srcdir, name):

    src = '%s/%s'%(srcdir, name)
    tgt = name
    print "reading", src

    s = open(src).read()

    lines = s.split('\n')

    _lines = []
    flag = False
    for line in lines:

        if line.startswith("% CUT HERE"):

            flag = not flag

        elif flag:

            _lines.append(line)

            if line.startswith("\includegraphics"):
                grab(srcdir, line)


    print "writing %d lines to %s" % (len(_lines), tgt)

    s = '\n'.join(_lines)

    f = open(tgt, 'w')
    print >>f, s


if __name__=="__main__":

    process(sys.argv[1], sys.argv[2])


    



