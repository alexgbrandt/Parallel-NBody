

proc runSim {testbin Nlist dt tEnd theta seed algList mpiargs} {

	set tList [];
	set E0list [];
	set Endlist [];

	set fname [string cat [file root [file tail $testbin]] "-" $dt "-" $tEnd "-" $theta ".txt"]
	set fp [open $fname "w"];
	puts $fp "#${mpiargs}"
	puts $fp "Alg,N,E0,Eend,T";
	close $fp;

	foreach N $Nlist algChoice $algList {

		set res [catch {exec mpirun {*}${mpiargs} $testbin $N $dt $tEnd $theta $seed $algChoice 2>@1} data];
		if {$res} {
			puts "$data"
			continue;
		}

		set E0 0;
		set Eend 0;
		set t 0;
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

		set fp [open $fname "a"];
		# foreach N $Nlist E0 $E0list Eend $Endlist T $tList {
			puts $fp "$algChoice,$N,$E0,$Eend,$t"
		# }

		close $fp;
	}


}


if {$::argc > 7} {
	lassign $::argv testbin Nlist dt tEnd theta seed algChoice mpiargs

	puts "NList: $Nlist";
	puts "dt: $dt";
	puts "tEnd: $tEnd";
	puts "theta: $theta";
	puts "seed: $seed";
	puts "algChoice: $algChoice";
	puts "mpiargs: $mpiargs";

	runSim $testbin $Nlist $dt $tEnd $theta $seed $algChoice $mpiargs

} else {
	puts "USAGE: tclsh RunMPISims.tcl <NBody-Bin> <N-Bodies-List> <delta-T> <end-T> <theta-MAC> <seed> <alg-choice-list> <mpi-args>"
}

