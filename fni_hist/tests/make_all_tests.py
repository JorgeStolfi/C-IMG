#! /usr/bin/python3
# Last edited on 2025-04-09 06:12:34 by stolfi

import os, sys, subprocess

def Er(msg):
  sys.stderr.write(msg)
  # ......................................................................
  
def bash(cmd):
  # The {cmd} must be a string with a command in {bash} syntax. 
  # Executes {cmd} with "/bin/bash".
  subprocess.run([ cmd ], shell=True, executable="/bin/bash")
  # ......................................................................

def run_command(cmd):
  # The {cmd} must be a list where the first elem is an executable file name
  # and the other elements are the command line arguments to be passed 
  # to its as the {argv} array. 
  result = subprocess.run(cmd, text = True)
  print(result.stderr)
  print(result.stdout)
  assert result.returncode == 0, f"** {cmd[0]} failed - returned status = {result}"
  return
  # ......................................................................

progDir = ".."
prog = "fni_hist.sh"

def main():
  # Last edited on 2025-02-28 18:05:02 by stolfi

  # testing the "-channel" option:
  run_test( "test-1", "0", 0, False, 0.100, None,  None, None,  None, None,  None ) 
  run_test( "test-1", "1", 1, False, 0.100, None,  None, None,  None, None,  None )
  run_test( "test-1", "2", 2, False, 0.100, None,  None, None,  None, None,  None )

  # Testing the "-range" and "-excludeRange" options:
  run_test( "test-3", "A", 1, False, 0.250, None,  None, None,   None, None,   None )
  run_test( "test-3", "B", 1, False, 0.250, None,  None, None,   -10.0, +1.0,  None )
  run_test( "test-3", "C", 1, False, 0.250, None,  -20.0, +3.7,  None, None,   None )

  # Testing the "-logScale" option:
  run_test( "test-4", "A", 0, False, 0.250, None,  None, None,   None, None,  None )
  run_test( "test-4", "B", 0, True,  0.250, None,  None, None,   None, None,  None )

  # Testing the "-binRound" option:
  run_test( "test-4", "C", 0, False, 0.250, None,  None, None,   None, None,   -1)
  run_test( "test-4", "D", 0, False, 0.250, None,  None, None,   None, None,   +1)

  # Testing the "-bins" option:
  run_test( "test-5", "A", 0, False, None, 200,  None, None,   None, None,  None )
  run_test( "test-6", "A", 0, False, None, 200,  None, None,   None, None,  None )

  return 0
  # ......................................................................
  
def run_test(name, tag, channel, logScale, step, bins, vMin, vMax, eMin, eMax, binRound):
  cmd = [ f"{progDir}/{prog}", "-channel", f"{channel:d}", ]
  if logScale:
    cmd.append( "-logScale" )
  if step != None:
    cmd += [ "-step", f"{step}" ]
  if bins != None:
    cmd += [ "-bins", f"{bins:d}" ]
  if vMin != None or vMax != None:
    assert vMin != None and vMax != None, "inconsistent {vMin,vMax}"
    cmd += [ "-range", f"{vMin}", f"{vMax}" ]
  if eMin != None or eMax != None:
    assert eMin != None and eMax != None, "inconsistent {eMin,eMax}"
    cmd += [ "-excludeRange", f"{eMin}", f"{eMax}" ]
  if binRound != None:
    cmd += [ "-binRound", f"{binRound:d}" ]
    
  cmd += [ "<", f"in/{name}.fni", ">", f"out/{name}-{tag}-hist.eps" ]
  
  cmd = " ".join(cmd)
  Er(cmd + "\n")
    
  bash(cmd)

  Er("\n")
  
  return
  # ......................................................................

main()
