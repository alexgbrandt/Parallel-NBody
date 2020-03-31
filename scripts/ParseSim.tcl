

proc runSim {testbin Nlist dt tEnd theta} {

	set tList [];
	set E0list [];
	set Endlist [];
	foreach N $Nlist {

		set data [exec $testbin $N $dt $tEnd $theta 2>@1];

		set lines [split $data "\n"];
		foreach line $lines {
			if {[regexp {E0[ ]*:[ ]*([-0-9.]+)} $line -> E0] != 0} {
				lappend E0list $E0;
			}
			if {[regexp {Eend[ ]*:[ ]*([-0-9.]+)} $line -> Eend] != 0} {
				lappend Endlist $Eend;
			}
			if {[regexp {Elapsed Time[ ]*:[ ]*([-0-9.]+)} $line -> t] != 0} {
				lappend tList $t;
			}
		}

	}

	set fname [string cat [file root [file tail $testbin]] "-" $dt "-" $tEnd "-" $theta ".txt"]
	set fp [open $fname "w"];
	puts $fp "N, E0, Eend, T";
	foreach N $Nlist E0 $E0list Eend $Endlist T $tList {
		puts $fp "$N, $E0, $Eend, $T"
	}

	close $fp;

}




if {$::argc > 4} {
	lassign $::argv testbin Nlist dt tEnd theta

	runSim $testbin $Nlist $dt $tEnd $theta;

} else {
	puts "USAGE: tclsh ParseSim.tcl <NBody-Bin> <N-Bodies> <delta-T> <end-T> <theat-MAC>"
}