#!/usr/bin/env perl
use strict;
use warnings;

# Screen proximity of ionizable side-chain atoms in a PDB structure.
# A short distance is a flag for inspection, not a protonation-state prediction.
my ($file, $acid_arg, $basic_arg) = @ARGV;
die "Usage: perl ionizable_pairs.pl <pdb_file> [acidic_cutoff_A] [basic_cutoff_A]\n"
    unless defined $file && @ARGV <= 3;

sub cutoff {
    my ($value, $label) = @_;
    unless (defined $value) {
        print "Cut-off for $label residues (A)? ";
        $value = <STDIN>;
    }
    die "Invalid $label cutoff: enter a positive number\n"
        unless defined($value) && $value =~ /^\s*(?:\d+(?:\.\d*)?|\.\d+)\s*$/ && $value > 0;
    return $value + 0;
}
my $acid_cut = cutoff($acid_arg, 'acidic');
my $basic_cut = cutoff($basic_arg, 'basic');

my %selected = (
    ASP => [qw(OD1 OD2)], GLU => [qw(OE1 OE2)],
    LYS => [qw(NZ)], ARG => [qw(NE NH1 NH2)],
    HIS => [qw(ND1 NE2)],
);
open my $fh, '<', $file or die "Cannot open '$file': $!\n";
my (%residues, @order);
while (my $line = <$fh>) {
    last if $line =~ /^ENDMDL/;
    next unless substr($line, 0, 6) eq 'ATOM  ' && length($line) >= 54;
    my $alt = substr($line, 16, 1);
    next unless $alt eq ' ' || $alt eq 'A';
    my $atom = substr($line, 12, 4); $atom =~ s/\s+//g;
    my $name = substr($line, 17, 3); $name =~ s/\s+//g;
    next unless exists $selected{$name};
    my $chain = substr($line, 21, 1);
    my $seq = substr($line, 22, 4); $seq =~ s/\s+//g;
    my $icode = substr($line, 26, 1);
    my $key = join '|', $chain, $seq, $icode, $name;
    unless (exists $residues{$key}) {
        $residues{$key} = { name => $name, chain => $chain,
            seq => $seq, icode => $icode, atoms => {} };
        push @order, $key;
    }
    my @xyz = map { 0 + substr($line, $_, 8) } (30, 38, 46);
    $residues{$key}{atoms}{$atom} //= \@xyz;
}
close $fh;

sub label {
    my ($r) = @_;
    my $chain = $r->{chain} eq ' ' ? '(blank)' : $r->{chain};
    my $ins = $r->{icode} eq ' ' ? '' : $r->{icode};
    return "$r->{name} $r->{seq}$ins (chain $chain)";
}
my (@acid, @basic, @his);
for my $key (@order) {
    my $r = $residues{$key};
    my @available = grep { exists $r->{atoms}{$_} } @{$selected{$r->{name}}};
    unless (@available) {
        warn "Skipping ", label($r), ": no ionizable side-chain atoms found\n";
        next;
    }
    my @missing = grep { !exists $r->{atoms}{$_} } @{$selected{$r->{name}}};
    warn "Partial atom set for ", label($r), ": missing @missing\n" if @missing;
    $r->{available} = \@available;
    if ($r->{name} =~ /^(?:ASP|GLU)$/) { push @acid, $r }
    elsif ($r->{name} eq 'HIS') { push @his, $r }
    else { push @basic, $r }
}

sub closest {
    my ($a, $b) = @_;
    my ($best, $atom_a, $atom_b);
    for my $aa (@{$a->{available}}) {
        for my $bb (@{$b->{available}}) {
            my $sum = 0;
            for my $k (0 .. 2) {
                $sum += ($a->{atoms}{$aa}[$k] - $b->{atoms}{$bb}[$k]) ** 2;
            }
            my $d = sqrt($sum);
            if (!defined($best) || $d < $best) {
                ($best, $atom_a, $atom_b) = ($d, $aa, $bb);
            }
        }
    }
    return ($best, $atom_a, $atom_b);
}
sub report {
    my ($title, $left, $right, $cut, $same) = @_;
    my @hits;
    for my $i (0 .. $#$left) {
        for my $j (($same ? $i + 1 : 0) .. $#$right) {
            my ($d, $aa, $bb) = closest($left->[$i], $right->[$j]);
            push @hits, [$d, $left->[$i], $right->[$j], $aa, $bb] if $d <= $cut;
        }
    }
    @hits = sort { $a->[0] <=> $b->[0] } @hits;
    printf "\n%s (cutoff %.1f A)\n", $title, $cut;
    if (!@hits) { print "  No pairs found.\n"; return }
    for my $h (@hits) {
        printf "  %-26s <-> %-26s %6.2f A  (%s-%s)\n",
            label($h->[1]), label($h->[2]), $h->[0], $h->[3], $h->[4];
    }
}

print "\nPDB: $file\n";
print "Distances are minimum separations between selected ionizable side-chain atoms, not pKa predictions.\n";
report('ACIDIC-ACIDIC (ASP/GLU)', \@acid, \@acid, $acid_cut, 1);
report('BASIC-BASIC (ARG/LYS)', \@basic, \@basic, $basic_cut, 1);
report('HISTIDINE-HISTIDINE (protonation-dependent)', \@his, \@his, $basic_cut, 1);
report('HISTIDINE vs ARG/LYS (protonation-dependent)', \@his, \@basic, $basic_cut, 0);
report('HISTIDINE vs ASP/GLU (protonation-dependent)', \@his, \@acid, $acid_cut, 0);
print "\nInterpret with local environment, solvent exposure and hydrogen bonding before pdb2gmx.\n";
