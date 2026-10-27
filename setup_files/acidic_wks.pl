#use strict;

print "What is the name of your pdb file?" ;
my $infile1 = <STDIN>;

my $outfile = "rascript";

# the pdb file with the experimental B-values
open(IN_ONE, "<$infile1") || die "Cannot open $infile1 $!";	

# the second input file 
# open(IN_TWO, "<$infile2") || die "Cannot open $infile2 $!";

# the output file, contains the results of the rotations
#open(PDB, ">$outfile") || die "Cannot open $outfile $!";


my @acidic_array;
my @first_tokens;
my @first_file;
my @first_temp;
my @glu_array1;
my @asp_array1;
my @lys_array1;
my @arg_array1;
my @his_array1;

my @data;

my $atom;
my $atom_num;
my $atom_type;
my $res_name;
my $chain;
my $res_num;
my $x_coord;
my $y_coord;
my $z_coord;
my $a;
my $b;
my $c;
my $chain;
my $loop_counter = 0;
my $res_counter = 0;
my $at_counter = 0;

my @chain_A;
my @chain_B;
my @chain_C;
my @chain_D;
my @chain_E;

while (<IN_ONE>)	{
	chomp ($_);
	@first_tokens = "";
	@first_temp = "";
	
    if ($_ =~ /^ATOM/)	{	
	$loop_counter++;
	@first_tokens = split(/\s+/, $_);	# split each line at space
	push(@first_temp, $first_tokens[0]);
	push(@first_temp, $first_tokens[1]);
	push(@first_temp, $first_tokens[2]);
	push(@first_temp, $first_tokens[3]);
	push(@first_temp, $first_tokens[4]);	# stored the number associated with each residue
	push(@first_temp, $first_tokens[5]);	# store x-value into @first_temp
	push(@first_temp, $first_tokens[6]);	# store y-value into @first_temp
	push(@first_temp, $first_tokens[7]);	# store z-value into @first_temp
	push(@first_temp, $first_tokens[8]);
	push(@first_temp, $first_tokens[9]);
	push(@first_temp, $first_tokens[10]);
	push @first_file, [@first_temp];	# store @static_temp into @static_coordinates array	
    }	
	
	
}

close(IN_ONE);

#START OF GLU

#print " checking for residue => $first_file[0][7]\n";
for ($i=0;$i<=$#first_file; $i++) {
	if (($first_file[$i][4] eq "GLU")) {
		if ($first_file[$i][3] eq "OE1") {
	   		push (@glu_array1, $first_file[$i]);
		}
	}
}


#print ("acidic residues\n");
# printing out each row	
#for my $m (0..$#glu_array1)	{
#	print "$m is:  @{$glu_array1[$m]}\n";
#}

my @distance_array1;
my @temp;
my @acid_array;
my @glu_array2;

#READING IN FIRST O
for (my $i=0; $i<=$#glu_array1; $i++)	{
        
	my $xx1 = $glu_array1[$i][7];
        my $yy1 = $glu_array1[$i][8];
        my $zz1 = $glu_array1[$i][9];

         for (my $j=0;$j<=$#first_file; $j++) {
	 if (($first_file[$j][4] eq "GLU")) {
		if (($first_file[$j][3] eq "OE2")){
	   		push (@glu_array2, $first_file[$j]);
		}
	}
}



#READING IN SECOND O AND WORKING OUT COORDINATES HALWAY BETWEEN TWO Os


for (my $j=$i; $j<=$#glu_array1; $j++)	{

                my $xx2 = $glu_array2[$j][7];
		my $yy2 = $glu_array2[$j][8];
		my $zz2 = $glu_array2[$j][9];
		if ($glu_array1[$i][6] eq $glu_array2[$j][6]) {
		
		my $newx = ($xx2+$xx1) / 2;
                my $newy = ($yy2+$yy1) / 2;
                my $newz = ($zz2+$zz1) / 2;		

		
		$temp[0] = $glu_array1[$i][4];
		$temp[1] = $glu_array1[$i][3];
		$temp[2] = $glu_array1[$i][6];
		$temp[3] = $glu_array2[$j][4];
                $temp[4] = $glu_array2[$j][3];
		$temp[5] = $glu_array2[$j][6];
		$temp[6] = $newx;
		$temp[7] = $newy;
		$temp[8] = $newz;
		push @acid_array, [@temp];
		}
	}
	}


	

#for my $m (0..$#acid_array)	{
#        print "new coordinates at $m are:  @{$acid_array[$m]}\n";

#}

#print "checking x coord => $acid_array[1][8]\n";

#START OF ASP

#print " checking for residue => $first_file[0][7]\n";
for ($i=0;$i<=$#first_file; $i++) {
	if (($first_file[$i][4] eq "ASP")) {
		if ($first_file[$i][3] eq "OD1") {
	   		push (@asp_array1, $first_file[$i]);
		}
	}
}


my @temp;
my @acid_array1;
my @asp_array2;

for (my $i=0; $i<=$#asp_array1; $i++)	{
        
	my $xx1 = $asp_array1[$i][7];
        my $yy1 = $asp_array1[$i][8];
        my $zz1 = $asp_array1[$i][9];

         for (my $j=0;$j<=$#first_file; $j++) {
	 if (($first_file[$j][4] eq "ASP")) {
		if (($first_file[$j][3] eq "OD2")){
	   		push (@asp_array2, $first_file[$j]);
		}
	}
}







for (my $j=$i; $j<=$#asp_array1; $j++)	{

                my $xx2 = $asp_array2[$j][7];
		my $yy2 = $asp_array2[$j][8];
		my $zz2 = $asp_array2[$j][9];
		if ($asp_array1[$i][6] eq $asp_array2[$j][6]) {
		
		my $newxx = ($xx2+$xx1) / 2;
                my $newyy = ($yy2+$yy1) / 2;
                my $newzz = ($zz2+$zz1) / 2;		

		
		$temp[0] = $asp_array1[$i][4];
		$temp[1] = $asp_array1[$i][3];
		$temp[2] = $asp_array1[$i][6];
		$temp[3] = $asp_array2[$j][4];
                $temp[4] = $asp_array2[$j][3];
		$temp[5] = $asp_array2[$j][6];
		$temp[6] = $newxx;
		$temp[7] = $newyy;
		$temp[8] = $newzz;
		push @acid_array1, [@temp];
		}
	}
	}


	

#for my $n (0..$#acid_array1)	{
#        print "new coordinates at $n are:  @{$acid_array1[$n]}\n";

#}

#print "checking x coord => $acid_array1[1][8]\n";


#START OF LYS

#print " checking for residue => $first_file[0][7]\n";
for ($i=0;$i<=$#first_file; $i++) {
	if (($first_file[$i][4] eq "LYS")) {
		if ($first_file[$i][3] eq "N") {
	   		push (@lys_array1, $first_file[$i]);
		}
	}
}
#print ("acidic residues\n");
# printing out each row	
#for my $m (0..$#glu_array1)	{
#	print "$m is:  @{$glu_array1[$m]}\n";
#}


my @temp;
my @acid_array2;
my @lys_array2;

for (my $i=0; $i<=$#lys_array1; $i++)	{
        
	my $xx2 = $lys_array1[$i][7];
        my $yy2 = $lys_array1[$i][8];
        my $zz2 = $lys_array1[$i][9];

                $temp[0] = $lys_array1[$i][4];
		$temp[1] = $lys_array1[$i][3];
		$temp[2] = $lys_array1[$i][6];
		$temp[3] = $xx2;
		$temp[4] = $yy2;
		$temp[5] = $zz2;
		push @lys_array2, [@temp];
         
		}





#for my $n (0..$#lys_array2)	{
#        print "new coordinates at $n are:  @{$lys_array2[$n]}\n";

#}

#print "checking x coord => $lys_array2[1][5]\n";


#START OF ARG

#print " checking for residue => $first_file[0][7]\n";
for ($i=0;$i<=$#first_file; $i++) {
	if (($first_file[$i][4] eq "ARG")) {
		if ($first_file[$i][3] eq "NH1") {
	   		push (@arg_array1, $first_file[$i]);
		}
	}
}
#print ("acidic residues\n");
# printing out each row	
#for my $m (0..$#glu_array1)	{
#	print "$m is:  @{$glu_array1[$m]}\n";
#}


my @temp;
my @basic_array;
my @arg_array2;

for (my $i=0; $i<=$#arg_array1; $i++)	{
        
	my $xx3 = $arg_array1[$i][7];
        my $yy3 = $arg_array1[$i][8];
        my $zz3 = $arg_array1[$i][9];

         for (my $j=0;$j<=$#first_file; $j++) {
	 if (($first_file[$j][4] eq "ARG")) {
		if (($first_file[$j][3] eq "NH2")){
	   		push (@arg_array2, $first_file[$j]);
		}
	}
}







for (my $j=$i; $j<=$#arg_array1; $j++)	{

                my $xx4 = $arg_array2[$j][7];
		my $yy4 = $arg_array2[$j][8];
		my $zz4 = $arg_array2[$j][9];
		if ($arg_array1[$i][6] eq $arg_array2[$j][6]) {
		
		my $newxxx = ($xx3+$xx4) / 2;
                my $newyyy = ($yy3+$yy4) / 2;
                my $newzzz = ($zz3+$zz4) / 2;		

		
		$temp[0] = $arg_array1[$i][4];
		$temp[1] = $arg_array1[$i][3];
		$temp[2] = $arg_array1[$i][6];
		$temp[3] = $arg_array2[$j][4];
                $temp[4] = $arg_array2[$j][3];
		$temp[5] = $arg_array2[$j][6];
		$temp[6] = $newxxx;
		$temp[7] = $newyyy;
		$temp[8] = $newzzz;
		push @basic_array, [@temp];
		}
	}
	}
	
	
	
	
	
#for my $n (0..$#basic_array)	{
#        print "new coordinates at $n are:  @{$basic_array[$n]}\n";

#}

#print "checking x coord => $basic_array[1][8]\n";


#START OF HIS



#print " checking for residue => $first_file[0][7]\n";
for ($i=0;$i<=$#first_file; $i++) {
	if (($first_file[$i][4] eq "HIS")) {
		if ($first_file[$i][3] eq "ND1") {
	   		push (@his_array1, $first_file[$i]);
		}
	}
}
#print ("acidic residues\n");
# printing out each row	
#for my $m (0..$#glu_array1)	{
#	print "$m is:  @{$glu_array1[$m]}\n";
#}


my @temp;

my @his_array2;
my @his_array3;

for (my $i=0; $i<=$#his_array1; $i++)	{
        
	my $xx5 = $his_array1[$i][7];
        my $yy5 = $his_array1[$i][8];
        my $zz5 = $his_array1[$i][9];

        	$temp[0] = $his_array1[$i][4];
		$temp[1] = $his_array1[$i][3];
		$temp[2] = $his_array1[$i][6];
		$temp[3] = $xx5;
		$temp[4] = $yy5;
		$temp[5] = $zz5;
		push @his_array2, [@temp];
		}
#	}
#}

for ($i=0;$i<=$#first_file; $i++) {
	if (($first_file[$i][4] eq "HIS")) {
		if ($first_file[$i][3] eq "NE2") {
	   		push (@his_array3, $first_file[$i]);
		}
	}
}
my @his_array4;

for (my $j=0; $j<=$#his_array3; $j++)	{
        
	my $xx6 = $his_array3[$j][7];
        my $yy6 = $his_array3[$j][8];
        my $zz6 = $his_array3[$j][9];

        	$temp[0] = $his_array3[$j][4];
		$temp[1] = $his_array3[$j][3];
		$temp[2] = $his_array3[$j][6];
		$temp[3] = $xx6;
		$temp[4] = $yy6;
		$temp[5] = $zz6;
		push @his_array4, [@temp];
		}

#for my $n (0..$#his_array2)	{
#        print "new coordinates at $n are:  @{$his_array2[$n]}\n";

#}
#}
#print "checking x coord => $his_array2[1][5]\n";

#for my $m (0..$#his_array4)	{
#        print "new coordinates at $m are:  @{$his_array4[$m]}\n";

#}

#CHECKING IF ELEMENTS OF ALL ARRAYS CAN BE ACCESSED

#print "checking x coord => $his_array4[1][3]\n";
#print "checking x coord => $his_array2[1][3]\n";
#print "checking x coord => $acid_array[1][3]\n";
#print "checking x coord => $acid_array1[1][3]\n";
#print "checking x coord => $lys_array2[1][3]\n";
#print "checking x coord => $basic_array[1][8]\n";

#COMBINING GLU and ASP to get one array: A_array

my @A_array;
my @temp;

for (my $i=0; $i<=$#acid_array; $i++)	{

                
		
		$temp[0] = $acid_array[$i][0];
		$temp[1] = "O";
		$temp[2] = $acid_array[$i][2];
		$temp[3] = $acid_array[$i][6];
                $temp[4] = $acid_array[$i][7];
		$temp[5] = $acid_array[$i][8];
#		$temp[6] = $acid_array[$i][6];
#		$temp[7] = $acid_array[$i][7];
#		$temp[8] = $acid_array[$i][8];

		
		push @A_array, [@temp];
	}




for (my $j=0; $j<=$#acid_array1; $j++)	{

                
		
		$temp[0] = $acid_array1[$j][0];
		$temp[1] = "O";
		$temp[2] = $acid_array1[$j][2];
		$temp[3] = $acid_array1[$j][6];
                $temp[4] = $acid_array1[$j][7];
		$temp[5] = $acid_array1[$j][8];
#		$temp[6] = $acid_array1[$j][6];
#		$temp[7] = $acid_array1[$j][7];
#		$temp[8] = $acid_array1[$j][8];

		
		push @A_array, [@temp];
	}


	

#for my $m (0..$#A_array)	{
#        print "new coordinates at $m are:  @{$A_array[$m]}\n";

#}
#CHECKING A_array
#print "checking x coord => $A_array[1][5]\n";
#print "checking x coord => $A_array[9][5]\n";

#COMBINING LYS, ARG and HIS to get one array: B_array

my @B_array;
my @temp;

for (my $i=0; $i<=$#lys_array2; $i++)	{

                
		
		$temp[0] = $lys_array2[$i][0];
		$temp[1] = $lys_array2[$i][1];
		$temp[2] = $lys_array2[$i][2];
		$temp[3] = $lys_array2[$i][3];
                $temp[4] = $lys_array2[$i][4];
		$temp[5] = $lys_array2[$i][5];
		

		
		push @B_array, [@temp];
	}
	


for (my $j=0; $j<=$#basic_array; $j++)	{

                
		
		$temp[0] = $basic_array[$j][0];
		$temp[1] = "N";
		$temp[2] = $basic_array[$j][2];
		$temp[3] = $basic_array[$j][6];
                $temp[4] = $basic_array[$j][7];
		$temp[5] = $basic_array[$j][8];
#		$temp[6] = $basic_array[$j][6];
#		$temp[7] = $basic_array[$j][7];
#		$temp[8] = $basic_array[$j][8];
		
		push @B_array, [@temp];
	}


	
for (my $k=0; $k<=$#his_array2; $k++)	{

                
		
		$temp[0] = $his_array2[$k][0];
		$temp[1] = $his_array2[$k][1];
		$temp[2] = $his_array2[$k][2];
		$temp[3] = $his_array2[$k][3];
                $temp[4] = $his_array2[$k][4];
		$temp[5] = $his_array2[$k][5];
		$temp[6] = $his_array2[$k][6];
		$temp[7] = $his_array2[$k][7];
		$temp[8] = $his_array2[$k][8];
		
		push @B_array, [@temp];
	}


for (my $l=0; $l<=$#his_array4; $l++)	{

                
		
		$temp[0] = $his_array4[$l][0];
		$temp[1] = $his_array4[$l][1];
		$temp[2] = $his_array4[$l][2];
		$temp[3] = $his_array4[$l][3];
                $temp[4] = $his_array4[$l][4];
		$temp[5] = $his_array4[$l][5];
		$temp[6] = $his_array4[$l][6];
		$temp[7] = $his_array4[$l][7];
		$temp[8] = $his_array4[$l][8];
		
		push @B_array, [@temp];
	}
	
#for my $m (0..$#B_array)	{
#        print "new coordinates at $m are:  @{$B_array[$m]}\n";
	
	
#}	

#CHECKING B_ARRAY
#print "checking x coord => $B_array[0][0]\n";
#print "checking x coord => $B_array[0][1]\n";
#print "checking x coord => $B_array[0][2]\n";
#print "checking x coord => $B_array[0][3]\n";
#print "checking x coord => $B_array[0][4]\n";
#print "checking x coord => $B_array[0][5]\n";
#print "checking x coord => $B_array[1][0]\n";
#print "checking x coord => $B_array[1][1]\n";
#print "checking x coord => $B_array[1][2]\n";
#print "checking x coord => $B_array[1][3]\n";
#print "checking x coord => $B_array[1][4]\n";
#print "checking x coord => $B_array[1][4]\n";

#DISTANCE CALCULATION FOR A_array

print "Cut-off for acidic residues?" ;
my $acut = <STDIN>;


my @distance_array;
my @temp;
for (my $i=0; $i<=$#A_array; $i++)	{

	my $x1 = $A_array[$i][3];
	my $y1 = $A_array[$i][4];
	my $z1 = $A_array[$i][5];
	
	for (my $j=$i+1; $j<=$#A_array; $j++)	{
	
		my $x2 = $A_array[$j][3];
		my $y2 = $A_array[$j][4];
		my $z2 = $A_array[$j][5];
		
		my $distance = sqrt((($x2-$x1)*($x2-$x1)) + (($y2-$y1)*($y2-$y1)) + (($z2-$z1)*($z2-$z1)));
		if ($distance<=$acut) { 
		$temp[0] = $A_array[$i][0];
		$temp[1] = $A_array[$i][1];
		$temp[2] = $A_array[$i][2];
		$temp[3] = $A_array[$j][0];
		$temp[4] = $A_array[$j][1];
                $temp[5] = $A_array[$j][2];
		$temp[6] = $distance;
		push @distance_array, [@temp];
		}
	}

}

for my $m (0..$#distance_array)	{
        print "acidic residues in close proximity $m are:  @{$distance_array[$m]}\n";
	
	
}

#DISTANCE CALCULATION FOR B_array

print "Cut-off for basic residues?" ;
my $bcut = <STDIN>;

my @distance_array2;
my @temp;
for (my $i=0; $i<=$#B_array; $i++)	{

	my $x1 = $B_array[$i][3];
	my $y1 = $B_array[$i][4];
	my $z1 = $B_array[$i][5];
	
	for (my $j=$i+1; $j<=$#B_array; $j++)	{
	
		my $x2 = $B_array[$j][3];
		my $y2 = $B_array[$j][4];
		my $z2 = $B_array[$j][5];
		
		my $distance = sqrt((($x2-$x1)*($x2-$x1)) + (($y2-$y1)*($y2-$y1)) + (($z2-$z1)*($z2-$z1)));
		if ($distance<=$bcut) { 
		if (($B_array[$i][2] ne $B_array[$j][2])) {
		
		$temp[0] = $B_array[$i][0];
		$temp[1] = $B_array[$i][1];
		$temp[2] = $B_array[$i][2];
		$temp[3] = $B_array[$j][0];
		$temp[4] = $B_array[$j][1];
                $temp[5] = $B_array[$j][2];
		$temp[6] = $distance;
		push @distance_array2, [@temp];
		}
	}

}
}

for my $m (0..$#distance_array2)	{
        print "basic residues in close proximity $m are:  @{$distance_array2[$m]}\n";
	
	

}

#print "checking x coord => $distance_array2[1][0]\n";
#print "checking x coord => $distance_array2[1][2]\n";
#open(PDB, ">$outfile") || die "Cannot open $outfile $!";
for my $n (0..$#distance_array)	{
for my $m (0..$#distance_array2)	{
#printf PDB "select $distance_array2[1][0] $distance_array2[2][1]\n";
}
}
#system ("/shared/glea/app/rasmol/rasmol_32BIT -script $outfile $infile1");
