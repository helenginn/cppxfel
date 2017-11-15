from cppxfel import *

def run(argc, argv):
  totalLine = ' '.join(argv)
  runCommandLineArgs(totalLine)

  return

if __name__ == '__main__':
  import sys

  import signal
  signal.signal(signal.SIGINT, signal.SIG_DFL)

  run(len(sys.argv), sys.argv[0:])
