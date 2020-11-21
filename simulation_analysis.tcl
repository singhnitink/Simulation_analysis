 ##Keep this file in the simulation directory
 ## start VMD
 ## Load the psf and dcd file
 ## source simulation_analysis.tcl
 ## Analysing: Waters, SASA, ROG, RMSD, H-bonds, RMSF
 ##******************
set nframes [molinfo top get numframes]
puts "*********************************"
puts "Total number of frames is ${nframes}"
puts "*********************************"
puts "Counting number of waters within 10A of protein"
set sel [atomselect top "protein"]
set waters [atomselect top "(resname TIP3 and name OH2) and (pbwithin 10 of protein)"]

 ##Fitting the protin molecule
set frame0 [atomselect top "protein and name CA" frame 0]
for {set ix 1 } {$ix < $nframes } { incr ix } {
		set selx [atomselect top "protein and name CA" frame $ix]
		set all [atomselect top all frame $ix]
		$all move [measure fit $selx $frame0]
	}
#**Specifying file names to store data	
set outfile1 [open ".water_data.csv" w]
set outfile2 [open "sasa.csv" w]
set outfile3 [open "radius_of_gyration.csv" w]
set outfile4 [open "rmsd.csv" w]
set outfile5 [open "hbonds.csv" w]
set outfile6 [open "rmsf.csv" w]
#**
##**Calculation of waters, sasa, rog,rmsd

for {set i 0} {$i < $nframes} {incr i} {
		$waters frame $i
		$waters update
		$sel frame $i
		$sel update
		set waternumbers [$waters num]
	    set sasa [measure sasa 1.4 $sel]
		set rog [measure rgyr $sel]
		set selrms [atomselect top "protein and backbone" frame $i]
		set rms [measure rmsd $selrms $frame0]
		puts $outfile1 [format "%d %2.2f" $i $waternumbers]
		puts $outfile2 [format "%d %2.2f" $i $sasa]
		puts $outfile3 [format "%d %2.2f" $i $rog]
		puts $outfile4 [format "%d %2.2f" $i $rms]
		
	}

##**Calculation of h-bonds	
set D [atomselect top "protein and name N"]
set A [atomselect top "protein and name O"]	
for {set j 0} {$j < $nframes} {incr j} {
		$D frame $j
		$A frame $j
		$A update
		$D update
		set hbondList [measure hbonds 3.0 20 $D $A]
		set hbondNumber [llength [lindex $hbondList 0]]
		puts $outfile5 [format "%d %d" $j $hbondNumber]
	}

##**Calculation of RMSF	
set sel1 [atomselect top "protein and name CA"]
 for {set k 0} {$k < [$sel1 num]} {incr k} {
 set rmsf [measure rmsf $sel1 ]
 puts $outfile6 "[expr {$k+1}] [lindex $rmsf $k]"
 }	
close $outfile1
close $outfile2
close $outfile3	
close $outfile4
close $outfile5
close $outfile6
}
cd ..
}
puts "This script calculated SASA, waters within 10 A of protein, ROG, RMSD, RMSF and h-bonds"
puts "GREAT!! IT IS ALL-DONE"
exit