#! /usr/bin/gawk -f
# Last edited on 2017-07-10 19:33:57 by stolfilocal

BEGIN { 
  debug = 0;  # Should it print debugging state?
  do_filter = 1; # Should it apply the filter?
  state = 0;  
    # If 0, just copying stream.
    # If 1, deleting irrelevant output.
    # If 2, collecting continuation lines of error message text {emess}. 
  if (debug != 0) { printf "@@ state = %d --starting--\n", state > "/dev/stderr"; }
}

(do_filter == 0){ 
  # Pass everything through:
  print; 
  next;
}

(debug != 0){
  printf "@@ state = %d %s\n", state, $0 > "/dev/stderr";
}

/^[!][!]/ { 
  # Debug message line:
  print; 
  next;
}

/^[ ]/ { 
  # Possible continuation line:
  if (state == 2) { emess = (emess $0); next; }
}

// { 
  # Not a continuation line. If gathering error message, dump it:
  if (state == 2) { 
    if ($0 ~ /^[ ]/) { prog_error(("boh?")); }
    printf "%s:%s: ** %s\n", efile, eline, emess > "/dev/stderr";
    state = 0;
  }
}
 
/^Persistence of Vision[(]tm[)]/ { state = 1; next; }
/^Output Options/ { state = 0; }
/^Total Scene Processing Times/ { state = 1; next; }

/^File[ ]+[']*[^' ]+['][ ]+line */ {
  # Preamble of context lines and error message:
  efile=$2; gsub(/[']/, "", efile);
  eline=$4; gsub(/[:]/, "", eline);
  if (! match($0, /[:]/)) { prog_error(("no colon")); }
  emess = substr($0, RSTART+1); gsub(/^[ ]+/, "", emess);
  state = 2; 
  next;
}

(state == 0) {
  # Generic random line:
  print > "/dev/stderr"; 
  next;
}

(state == 1) { 
  # Undesired noise line:
  next;
}

(state == 2) {
  # Should have been handled before:
  prog_error(("boh?")); 
}
 

function prog_error(msg){
  printf "!! povray-output-filter: bug at output line %d - %s\n", FNR, msg > "/dev/stderr";
  print > "/dev/stderr";
  next;
}
