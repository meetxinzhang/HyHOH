#!/usr/bin/perl -w

use strict;
use warnings;
use Cwd;


my $dir = getcwd.'/pull_out';  # get path

opendir my $dh, $dir  or die "Can't open $dir: $!";
my $num_files = grep { -f "$dir/$_" } readdir($dh);

rewinddir($dh);  # so that it can read the dir again
my $num_gro = grep { /.gro/ } readdir($dh);  # chose .gro

my $n_conf = $num_gro;
# print $num_files;
# print $dir;



# loop g_dist command - measure distance in each frame, write to a file
for (my $i=0; $i<=$n_conf; $i++) {
    print "Processing configuration $i...\n";
    system("gmx distance -s pull.tpr -f pull_out/conf${i}.gro -n index.ndx -oav dist${i}.xvg -select \'com of group \"receptor\" plus com of group \"ligand\"\' &>/dev/null");
}

# write output to single file
open(OUT, ">>summary_distances.dat");

for (my $j=0; $j<=$n_conf; $j++) {
    open(IN, "<dist${j}.xvg");
    my @array = <IN>;

    my $distance;
    # my $x;
    # my $y;
    # my $z;

    foreach $_ (@array) {
        if ($_ =~ /[#@]/) {
            # do nothing, it's a comment or formatting line
        } else {
            my @line = split(" ", $_);
            # $x = $line[1];
            # $y = $line[2];
            # $z = $line[3];
            # $distance = sqrt($x^2 + $y^2 + $z^2)
            $distance = $line[1]
            # Code above(2nd norm calculating) does the same thing as -oav option
        }
    }

    close(IN);
    print OUT "$j\t$distance\n";
}

close(OUT);

# clean up
print "Cleaning up...\n";

for (my $k=0; $k<=$n_conf; $k++) {
    unlink "dist${k}.xvg";
}

exit;
