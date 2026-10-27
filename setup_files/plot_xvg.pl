#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename qw(fileparse);

@ARGV == 1 or die "Usage: perl plot_xvg.pl filename.xvg\n";
my $input = $ARGV[0];
open my $in, '<', $input or die "Cannot open $input: $!\n";
my ($title, $xlabel, $ylabel) = ('GROMACS data', 'Time (ps)', 'Value');
my $legend;
my (@x, @y);
while (my $line = <$in>) {
    if ($line =~ /^\@\s+title\s+"(.*?)"/)  { $title = $1; next; }
    if ($line =~ /^\@\s+xaxis\s+label\s+"(.*?)"/) { $xlabel = $1; next; }
    if ($line =~ /^\@\s+s0\s+legend\s+"(.*?)"/) { $legend = $1; next; }
    next if $line =~ /^\s*[\@#&]/ || $line =~ /^\s*$/;
    my @v = split ' ', $line;
    next unless @v >= 2 && $v[0] =~ /^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?$/ && $v[1] =~ /^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?$/;
    push @x, 0 + $v[0]; push @y, 0 + $v[1];
}
close $in;
# GROMACS energy XVG files may carry a generic, incorrect y-axis unit.
# Prefer the first series legend and deliberately omit units.
if (defined $legend && length $legend) { $title = $legend; $ylabel = $legend; }
else { $ylabel = $title eq "GROMACS Energies" ? "Value" : $title; }
@x >= 2 or die "Need at least two numeric data rows in $input\n";
my ($xmin, $xmax, $ymin, $ymax) = ($x[0], $x[0], $y[0], $y[0]);
for my $i (0 .. $#x) {
    $xmin = $x[$i] if $x[$i] < $xmin;
    $xmax = $x[$i] if $x[$i] > $xmax;
    $ymin = $y[$i] if $y[$i] < $ymin;
    $ymax = $y[$i] if $y[$i] > $ymax;
}
my $xr = $xmax - $xmin || 1;
my $yr = $ymax - $ymin || 1;
$ymin -= $yr * .06; $ymax += $yr * .06;
$yr = $ymax - $ymin;
my ($w, $h, $left, $top, $pw, $ph) = (960, 560, 105, 70, 800, 400);
sub esc { my $s = shift; $s =~ s/&/&amp;/g; $s =~ s/</&lt;/g; $s =~ s/>/&gt;/g; $s =~ s/"/&quot;/g; return $s; }
sub sx { return $left + ($_[0] - $xmin) / $xr * $pw; }
sub sy { return $top + $ph - ($_[0] - $ymin) / $yr * $ph; }
my ($stem, $dir) = fileparse($input, qr/\.xvg\z/i);
my $output = $dir . $stem . '.html';
open my $out, '>', $output or die "Cannot write $output: $!\n";
print $out '<!doctype html><html><head><meta charset="utf-8"><title>', esc($title), '</title><style>body{font-family:Arial,sans-serif;background:#fafafa;text-align:center}svg{max-width:100%;height:auto;background:white;border:1px solid #ddd}</style></head><body><svg viewBox="0 0 960 560" xmlns="http://www.w3.org/2000/svg">', "\n";
print $out '<text x="480" y="35" text-anchor="middle" font-size="22" font-weight="bold">', esc($title), "</text>\n";
for my $i (0..5) {
    my $xx = $left + $pw * $i / 5;
    my $yy = $top + $ph * $i / 5;
    printf $out '<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#eee"/><text x="%.1f" y="%d" font-size="13" text-anchor="middle">%.4g</text>', $xx,$top,$xx,$top+$ph,$xx,$top+$ph+23,$xmin+$xr*$i/5;
    printf $out '<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#eee"/><text x="%d" y="%.1f" font-size="13" text-anchor="end">%.4g</text>', $left,$yy,$left+$pw,$yy,$left-9,$yy+4,$ymax-$yr*$i/5;
}
print $out qq{<rect x="$left" y="$top" width="$pw" height="$ph" fill="none" stroke="#444"/>\n};
print $out '<polyline fill="none" stroke="#1768ac" stroke-width="1.5" stroke-linejoin="round" points="';
for my $i (0..$#x) { printf $out '%.2f,%.2f ', sx($x[$i]), sy($y[$i]); }
print $out '"/>', "\n";
print $out '<text x="505" y="535" text-anchor="middle" font-size="17">', esc($xlabel), "</text>\n";
print $out '<text transform="translate(24 270) rotate(-90)" text-anchor="middle" font-size="17">', esc($ylabel), "</text>\n";
print $out '</svg></body></html>', "\n";
close $out;
print "Created $output\n";
