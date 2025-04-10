#!/usr/bin/perl
use strict;
use warnings;

### add barocde (from R1 and R3) to .pairsam based on field
### format:
# CTACTAAGTTTAATCG:VH00454:1:AAC3FTHM5:1:1101:28835:1000  !       0       !       0       -       -       MM      CTACTAAGTTTAATCG:VH00454:1:AAC3FTHM5:1:1101:28835:100081chr2986668493131M=98666849-31TAAAAAACGTGAAAAATGTGAAACGCACACNCCCCCCCCCCCCCCCCCCCCCCC;CCCCCC#NM:i:2MD:Z:23T6T0MC:Z:31MAS:i:25XS:i:0CR:Z:CTACTAAGTTTAATCGYt:Z:MM    CTACTAAGTTTAATCG:VH00454:1:AAC3FTHM5:1:1101:28835:1000161chr2986668493431M=9866684931TAAAAAACGTGAAAAATGTGAAACGCACACTCCC-CCCCCCCCCCCCCCCCCCCCCCCC;CCNM:i:1MD:Z:23T7MC:Z:31MAS:i:26XS:i:19CR:Z:CTACTAAGTTTAATCGYt:Z:MM

open IN, "cat $ARGV[0]|" or die $!; ### dont print out ^Y
open OUT, " > $ARGV[0]\_bc";
my $bc1;
my $bc2;
my $corrupt = 0;

while(<IN>){
	chomp;
	if($_ =~ m/^#columns/){
		print OUT $_." barcode1 barcode2\n";
	}
	elsif($_ =~ m/^\#/){
		print OUT $_."\n";
	}
	else{
    	my @sp = split/\s+/, $_;
    	### field check: 
    	### 1. there should always be 10 fields before adding barcodes.
    	### 2. XX pairs need to be removed. Report a total of XX pairs removed.
    	### 3. Barcodes type is correct?
    	if($sp[7] eq "XX"){
    		$corrupt++;
    		next;  
    	}elsif(scalar(@sp) != 10){
    		$corrupt++;
    		next;
    	}
    	my @tmp1 = split /[:^]+/, $sp[8];
		$bc1 = $tmp1[0];
		my @tmp2 = split /[:^]+/, $sp[9];
		$bc2 = $tmp2[0];
    	my $output = $sp[0];
    	if(length($bc1) != 16 || length($bc2) != 16){
    		die "Barcodes at head of readname is not 16bp. Are you sure you are using the right pairsam?\n";
    	}
    	foreach my $i (1 .. $#sp){
    		$output  .= "\t$sp[$i]";
    	}
    	$output  .= "\t$bc1";
    	$output  .= "\t$bc2";
    	print OUT $output."\n";
    	}
}
close IN;
close OUT;

print "Pairs removed due to unmapped (type: XX): ", $corrupt, "\n";
system("mv $ARGV[0]\_bc $ARGV[0]");






