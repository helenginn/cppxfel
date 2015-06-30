from cppxfel import *

def run(argc, argv):
  cppxfelMain(argc, argv)

  return

if __name__ == '__main__':
  import sys
  run(len(sys.argv), sys.argv[0:])
