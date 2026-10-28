#!/usr/bin/env perl
use strict;
use warnings;

# Preserve the screening geometry of the original acidic_wks.pl.
# These are proximity flags, not necessarily side-chain charge contacts.
my ($file, $acid_cut, $basic_cut) = @ARGV;
die "Usage: perl ionizable_pairs.pl <pdb_file> [acidic_cutoff_A] [basic_cutoff_A]\n"
    unless defined $file && @ARGV <= 3;

sub cutoff {
    my ($value, $label) = @_;
    if (!defined $value) {
        print "Cut-off for $label residues? ";
        $value = <STDIN>;
    }
    defined($value) && $value =~ /^\s*(?:\d+(?:\.\d*)?|\.\d+)\s*$/ && $value > 0
        or die "Invalid $label cutoff: expected a positive number in Angstroms\n";
    return $value + 0;
}

open my $fh, '<', $file or die "Cannot open $file: $!\n";
my (%res, @order);
while (my $line = <$fh>) {
    last if $line =~ /^ENDMDL/;
    next unless $line =~ /^ATOM  /;
    my $alt = substr($line, 16, 1);
    next unless $alt eq ' ' || $alt eq 'A';
    my $atom = substr($line, 12, 4); $atom =~ s/\s+//g;
    my $name = substr($line, 17, 3); $name =~ s/\s+//g;
    next unless $name =~ /^(?:ASP|GLU|LYS|ARG|HIS)$/;
    my $chain = substr($line, 21, 1);
    my $num = substr($line, 22, 4); $num =~ s/\s+//g;
    my $icode = substr($line, 26, 1);
    my $key = join ':', $chain, $num, $icode;
    if (!exists $res{$key}) {
        $res{$key} = { name => $name, chain => $chain, num => $num,
                        icode => $icode, atoms => {} };
        push @order, $key;
    }
    my @xyz = map { 0 + substr($line, $_, 8) } (30, 38, 46);
    $res{$key}{atoms}{$atom} //= \@xyz;
}
close $fh;

sub midpoint {
    my ($r, $a, $b) = @_;
    return unless exists $r->{atoms}{$a} && exists $r->{atoms}{$b};
    return [ map { ($r->{atoms}{$a}[$_] + $r->{atoms}{$b}[$_]) / 2 } 0..2 ];
}
sub point {
    my ($r, $atom) = @_;
    return $r->{atoms}{$atom};
}
sub distance {
    my ($a, $b) = @_;
    return sqrt(($a->[0]-$b->[0])**2 + ($a->[1]-$b->[1])**2 + ($a->[2]-$b->[2])**2);
}
sub label {
    my ($r, $atom) = @_;
    my $chain = $r->{chain} =~ /\S/ ? $r->{chain} : '(blank)';
    return sprintf('%s %s %s (chain %s)', $r->{name}, $r->{num}.$r->{icode}, $atom, $chain);
}

# Match original ordering: GLU then ASP; LYS then ARG then HIS ND1 then HIS NE2.
my (@acid, @basic);
for my $type (qw(GLU ASP)) {
    for my $key (@order) {
        my $r = $res{$key}; next unless $r->{name} eq $type;
        my $p = $type eq 'GLU' ? midpoint($r, 'OE1', 'OE2') : midpoint($r, 'OD1', 'OD2');
        push @acid, [$r, 'O', $p] if $p;
    }
}
for my $type (qw(LYS ARG)) {
    for my $key (@order) {
        my $r = $res{$key}; next unless $r->{name} eq $type;
        my $p = $type eq 'LYS' ? point($r, 'N') : midpoint($r, 'NH1', 'NH2');
        push @basic, [$r, 'N', $p] if $p;
    }
}
for my $atom (qw(ND1 NE2)) {
    for my $key (@order) {
        my $r = $res{$key}; next unless $r->{name} eq 'HIS';
        my $p = point($r, $atom);
        push @basic, [$r, $atom, $p] if $p;
    }
}

$acid_cut = cutoff($acid_cut, 'acidic');
$basic_cut = cutoff($basic_cut, 'basic');
print "PDB: $file\n";
print "Distances use original-script reference atoms (LYS backbone N); proximity is not proof of ionizable side-chain interaction.\n\n";

sub report_pairs {
    my ($title, $items, $cutoff, $skip_same) = @_;
    printf "%s (cutoff %.1f A)\n", $title, $cutoff;
    my $count = 0;
    for my $i (0..$#$items) {
        for my $j ($i+1..$#$items) {
            my ($r1, $a1, $p1) = @{$items->[$i]};
            my ($r2, $a2, $p2) = @{$items->[$j]};
            # Avoid contacts within the same residue, even for HIS ND1/NE2.
            next if $skip_same && $r1 == $r2;
            my $d = distance($p1, $p2);
            next if $d > $cutoff;
            printf "  %-30s <-> %-30s %6.2f A\n", label($r1,$a1), label($r2,$a2), $d;
            $count++;
        }
    }
    print "  No pairs found.\n" unless $count;
    print "\n";
}
report_pairs('ACIDIC RESIDUE PROXIMITY (ASP/GLU)', \@acid, $acid_cut, 1);
report_pairs('BASIC RESIDUE PROXIMITY (ARG/LYS/HIS)', \@basic, $basic_cut, 1);
print "Note: HIS is protonation-dependent. LYS N is a backbone atom, not NZ.\n";
