#!/usr/bin/env perl
use strict;
use warnings;

# Screen nearby ionizable side chains before choosing pdb2gmx protonation states.
# Distances are between representative side-chain charge centres, in Angstroms.
my $file = shift @ARGV or die "Usage: perl acidic_wks.pl protein.pdb [acidic_cutoff_A] [basic_cutoff_A]\n";
die "Too many arguments\n" if @ARGV > 2;

sub cutoff {
    my ($value, $prompt) = @_;
    unless (defined $value) {
        print "$prompt (Angstroms): ";
        $value = <STDIN>;
        die "No cutoff supplied\n" unless defined $value;
        chomp $value;
    }
    die "Cutoff must be a positive number\n"
        unless $value =~ /\A(?:\d+(?:\.\d*)?|\.\d+)\z/ && $value > 0;
    return $value;
}

my $acid_cut = cutoff(shift @ARGV, 'Cut-off for acidic residues');
my $basic_cut = cutoff(shift @ARGV, 'Cut-off for basic residues');

open my $fh, '<', $file or die "Cannot open '$file': $!\n";
my (%residues, @order);
while (my $line = <$fh>) {
    last if $line =~ /^ENDMDL\s*$/; # first model only
    next unless substr($line, 0, 6) eq 'ATOM  ';
    next if length($line) < 54;
    my $atom = substr($line, 12, 4); $atom =~ s/\s+//g;
    my $alt  = substr($line, 16, 1);
    next unless $alt eq ' ' || $alt eq 'A';
    my $res  = substr($line, 17, 3); $res =~ s/\s+//g;
    next unless $res =~ /\A(?:ASP|GLU|LYS|ARG|HIS)\z/;
    my $chain = substr($line, 21, 1);
    my $seq   = substr($line, 22, 4); $seq =~ s/\s+//g;
    my $icode = substr($line, 26, 1);
    my $key = join '|', $chain, $seq, $icode, $res;
    my @xyz = map { 0 + substr($line, $_, 8) } (30, 38, 46);
    if (!exists $residues{$key}) {
        $residues{$key} = { res => $res, chain => $chain, seq => $seq,
                             icode => $icode, atoms => {} };
        push @order, $key;
    }
    $residues{$key}{atoms}{$atom} //= \@xyz;
}
close $fh;

my %centres = (
    ASP => [qw(OD1 OD2)], GLU => [qw(OE1 OE2)],
    LYS => [qw(NZ)], ARG => [qw(NH1 NH2)],
    HIS => [qw(ND1 NE2)],
);
my (@acid, @basic, @his);
for my $key (@order) {
    my $r = $residues{$key};
    my @names = @{$centres{$r->{res}}};
    my @coords;
    for my $name (@names) {
        if (!exists $r->{atoms}{$name}) {
            warn "Skipping ", describe($r), ": missing $name\n";
            @coords = ();
            last;
        }
        push @coords, $r->{atoms}{$name};
    }
    next unless @coords == @names;
    my @centre = map {
        my $dim = $_;
        my $sum = 0;
        $sum += $_->[$dim] for @coords;
        $sum / @coords;
    } (0 .. 2);
    $r->{centre} = \@centre;
    if ($r->{res} =~ /\A(?:ASP|GLU)\z/) { push @acid, $r }
    elsif ($r->{res} eq 'HIS') { push @his, $r }
    else { push @basic, $r }
}

sub describe {
    my ($r) = @_;
    my $chain = $r->{chain} eq ' ' ? '(blank)' : $r->{chain};
    my $ins = $r->{icode} eq ' ' ? '' : $r->{icode};
    return "$r->{res} $r->{seq}$ins (chain $chain)";
}
sub distance {
    my ($a, $b) = @_;
    my $sum = 0;
    for my $k (0 .. 2) {
        $sum += ($a->{centre}[$k] - $b->{centre}[$k]) ** 2;
    }
    return sqrt($sum);
}
sub report {
    my ($title, $left, $right, $cut, $same_list) = @_;
    my @hits;
    for my $i (0 .. $#$left) {
        my $start = $same_list ? $i + 1 : 0;
        for my $j ($start .. $#$right) {
            my $d = distance($left->[$i], $right->[$j]);
            push @hits, [$d, $left->[$i], $right->[$j]] if $d <= $cut;
        }
    }
    @hits = sort { $a->[0] <=> $b->[0] } @hits;
    printf "\n%s (cutoff %.1f A)\n", $title, $cut;
    if (!@hits) { print "  No pairs found.\n"; return }
    for my $hit (@hits) {
        printf "  %-26s <-> %-26s %6.2f A\n",
            describe($hit->[1]), describe($hit->[2]), $hit->[0];
    }
}

print "\nPDB: $file\n";
print "Distances are approximate side-chain charge-centre separations, not pKa predictions.\n";
report('ACIDIC-ACIDIC (ASP/GLU)', \@acid, \@acid, $acid_cut, 1);
report('BASIC-BASIC (ARG/LYS)', \@basic, \@basic, $basic_cut, 1);
report('HISTIDINE-HISTIDINE (protonation-dependent)', \@his, \@his, $basic_cut, 1);
report('HISTIDINE vs ARG/LYS (protonation-dependent)', \@his, \@basic, $basic_cut, 0);
report('HISTIDINE vs ASP/GLU (protonation-dependent)', \@his, \@acid, $acid_cut, 0);
print "\nInterpret with local environment, solvent exposure and hydrogen bonding before pdb2gmx.\n";
