#!/usr/bin/env perl
use strict;
use warnings;
use Time::HiRes qw(time);
use POSIX qw(strftime);
use Cwd qw(abs_path getcwd);
use File::Basename qw(dirname);
use File::Path qw(make_path);

my $manifest = shift @ARGV // "bench/bam-shape/inputs.tsv";
my $outdir = $ENV{OUTDIR} // "bench/bam-shape/results/htslib-test-view";
my $repeats = $ENV{REPEATS} // 5;
my $loops = $ENV{LOOPS} // 1;
my @threads = split /\s+/, ($ENV{THREADS_LIST} // "1 2 4 8");
my @modes = split /\s+/, ($ENV{MODES} // "decode readwrite");
my $base = $ENV{BASE_TEST_VIEW} // "../../htslib/test/test_view";
my $product = $ENV{PRODUCT_TEST_VIEW} // "../htslib-bam-product/test/test_view";
my $samtools = $ENV{SAMTOOLS} // "../../samtools/samtools";

make_path($outdir);

open my $tfh, ">", "$outdir/timings.tsv" or die "timings.tsv: $!\n";
open my $mfh, ">", "$outdir/metadata.tsv" or die "metadata.tsv: $!\n";
open my $cfh, ">", "$outdir/checks.tsv" or die "checks.tsv: $!\n";

print $tfh join("\t", qw(timestamp exe mode threads input shape region repeat loops seconds status command)), "\n";
print $mfh join("\t", qw(timestamp key value)), "\n";
print $cfh join("\t", qw(input shape region count status)), "\n";

sub ts {
    return strftime("%Y-%m-%dT%H:%M:%SZ", gmtime());
}

sub meta {
    my ($key, $value) = @_;
    $value = "" unless defined $value;
    $value =~ s/\t/ /g;
    print $mfh join("\t", ts(), $key, $value), "\n";
}

sub run_capture {
    my (@cmd) = @_;
    my $out = qx{@cmd 2>&1};
    my $status = $? == 0 ? "ok" : "fail";
    chomp $out;
    return ($status, $out);
}

sub cmd_string {
    return join " ", map {
        my $s = $_;
        $s =~ s/'/'\\''/g;
        $s =~ /[^A-Za-z0-9_.,:\/=+-]/ ? "'$s'" : $s;
    } @_;
}

sub test_view_cmd {
    my ($exe, $mode, $threads, $path, $region) = @_;
    my @cmd = ($exe, "-@", $threads);
    if ($mode eq "decode") {
        push @cmd, "-B", "-p", "/dev/null";
    } elsif ($mode eq "readwrite") {
        push @cmd, "-b", "-p", "/dev/null";
    } else {
        die "unknown mode: $mode\n";
    }
    push @cmd, $path;
    push @cmd, $region if defined($region) && $region ne "-";
    return @cmd;
}

sub record_binary {
    my ($label, $path) = @_;
    meta("$label.path", $path);
    meta("$label.abspath", abs_path($path) // "missing");
    if (-x $path) {
        my ($status, $out) = run_capture($path);
        my ($first) = grep { /\S/ } split /\n/, $out;
        meta("$label.version_probe_status", $status);
        meta("$label.version_probe_first_line", $first // "");
    } else {
        meta("$label.version_probe_status", "missing");
    }
}

sub git_meta {
    my ($label, $repo) = @_;
    if (-d $repo) {
        my ($st1, $commit) = run_capture("git", "-C", $repo, "rev-parse", "HEAD");
        my ($st2, $status) = run_capture("git", "-C", $repo, "status", "--short");
        meta("$label.git_commit_status", $st1);
        meta("$label.git_commit", $commit);
        meta("$label.git_status_status", $st2);
        $status =~ s/\n/; /g;
        meta("$label.git_status_short", $status);
    }
}

meta("cwd", getcwd());
meta("manifest", $manifest);
meta("outdir", $outdir);
meta("repeats", $repeats);
meta("loops", $loops);
meta("threads_list", join(" ", @threads));
meta("modes", join(" ", @modes));
record_binary("base_test_view", $base);
record_binary("product_test_view", $product);
record_binary("samtools", $samtools);
git_meta("base_htslib", dirname(dirname($base)));
git_meta("product_htslib", dirname(dirname($product)));
my ($ld_status, $ld_version) = run_capture("pkg-config", "--modversion", "libdeflate");
meta("libdeflate.pkg_config_status", $ld_status);
meta("libdeflate.pkg_config_version", $ld_version);

open my $ifh, "<", $manifest or die "$manifest: $!\n";
my $header = <$ifh>;
while (my $line = <$ifh>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    my ($name, $path, $region, $shape) = split /\t/, $line;
    next unless defined $name && length $name;
    unless (-r $path) {
        print STDERR "missing input: $name ($path)\n";
        next;
    }

    my @count_cmd = ($samtools, "view", "-c", $path);
    push @count_cmd, $region if defined($region) && $region ne "-";
    my ($count_status, $count_out) = run_capture(@count_cmd);
    my ($count) = $count_out =~ /^(\d+)$/m;
    print $cfh join("\t", $name, $shape, $region, $count // "NA", $count_status), "\n";

    for my $threads (@threads) {
        for my $mode (@modes) {
            for my $pair ([base => $base], [product => $product]) {
                my ($label, $exe) = @$pair;
                next unless -x $exe;
                my @cmd = test_view_cmd($exe, $mode, $threads, $path, $region);
                my $cmd_string = cmd_string(@cmd);
                for my $repeat (1 .. $repeats) {
                    my $start = time();
                    my $status = "ok";
                    for (1 .. $loops) {
                        system @cmd;
                        if ($? != 0) {
                            $status = "fail";
                            last;
                        }
                    }
                    my $seconds = time() - $start;
                    print $tfh join("\t", ts(), $label, $mode, $threads, $name,
                                    $shape, $region, $repeat, $loops,
                                    sprintf("%.6f", $seconds), $status, $cmd_string), "\n";
                }
            }
        }
    }
}

close $ifh;
close $tfh;
close $mfh;
close $cfh;
print "$outdir\n";
