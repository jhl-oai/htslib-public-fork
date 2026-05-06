#!/usr/bin/env perl
use strict;
use warnings;
use Time::HiRes qw(time);
use POSIX qw(strftime);
use Cwd qw(abs_path getcwd);
use File::Path qw(make_path);
use IPC::Open3;
use Symbol qw(gensym);

my $manifest = shift @ARGV // "bench/bam-shape/inputs.tsv";
my $outdir = $ENV{OUTDIR} // "bench/bam-shape/results/bgzf-qsize-matrix";
my $exe = $ENV{TEST_VIEW} // "../htslib-bam-product/test/test_view";
my $samtools = $ENV{SAMTOOLS} // "../../samtools/samtools";
my $repeats = $ENV{REPEATS} // 5;
my $loops = $ENV{LOOPS} // 1;
my $warmups = $ENV{WARMUPS} // 1;
my @threads = split /\s+/, ($ENV{THREADS_LIST} // "4 8");
my @modes = split /\s+/, ($ENV{MODES} // "decode readwrite");
my @qsizes = split /\s+/, ($ENV{QSIZE_LIST} // "default 32 64");
my %input_filter = map { $_ => 1 } grep { length } split /\s+/, ($ENV{INPUTS} // "");

make_path($outdir);

open my $tfh, ">", "$outdir/timings.tsv" or die "timings.tsv: $!\n";
open my $mfh, ">", "$outdir/metadata.tsv" or die "metadata.tsv: $!\n";
open my $cfh, ">", "$outdir/checks.tsv" or die "checks.tsv: $!\n";

print $tfh join("\t", qw(timestamp mode threads qsize input shape region repeat loops seconds max_rss_bytes status command)), "\n";
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
    my ($mode, $threads, $path, $region) = @_;
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

sub run_timed_once {
    my ($qsize, @cmd) = @_;
    my $err = gensym;
    my $out = gensym;
    my $start = time();
    local $ENV{HTS_BGZF_QUEUE_SIZE};
    if ($qsize eq "default") {
        delete $ENV{HTS_BGZF_QUEUE_SIZE};
    } else {
        $ENV{HTS_BGZF_QUEUE_SIZE} = $qsize;
    }
    my $pid = open3(undef, $out, $err, "/usr/bin/time", "-l", @cmd);
    my $stderr = do { local $/; <$err> };
    my $stdout = do { local $/; <$out> };
    waitpid($pid, 0);
    my $seconds = time() - $start;
    my $status = $? == 0 ? "ok" : "fail";
    my ($rss) = $stderr =~ /(\d+)\s+maximum resident set size/;
    if ($status ne "ok" && defined($stdout) && length($stdout)) {
        $stderr .= "\nstdout:\n$stdout";
    }
    return ($status, $seconds, $rss // "NA", $stderr);
}

sub run_timed_loop {
    my ($qsize, @cmd) = @_;
    my $total = 0.0;
    my $max_rss = 0;
    my $status = "ok";
    my $diag = "";
    for (1 .. $loops) {
        my ($one_status, $seconds, $rss, $stderr) = run_timed_once($qsize, @cmd);
        $total += $seconds;
        $max_rss = $rss if $rss ne "NA" && $rss > $max_rss;
        if ($one_status ne "ok") {
            $status = "fail";
            $diag = $stderr;
            last;
        }
    }
    return ($status, $total, $max_rss || "NA", $diag);
}

meta("cwd", getcwd());
meta("manifest", $manifest);
meta("outdir", $outdir);
meta("test_view.path", $exe);
meta("test_view.abspath", abs_path($exe) // "missing");
meta("samtools.path", $samtools);
meta("samtools.abspath", abs_path($samtools) // "missing");
meta("repeats", $repeats);
meta("loops", $loops);
meta("warmups", $warmups);
meta("threads_list", join(" ", @threads));
meta("modes", join(" ", @modes));
meta("qsize_list", join(" ", @qsizes));
meta("inputs_filter", join(" ", sort keys %input_filter));

my ($git_status, $git_commit) = run_capture("git", "-C", "../htslib-bam-product", "rev-parse", "HEAD");
meta("product.git_commit_status", $git_status);
meta("product.git_commit", $git_commit);
my ($status_status, $git_short) = run_capture("git", "-C", "../htslib-bam-product", "status", "--short");
$git_short =~ s/\n/; /g;
meta("product.git_status_status", $status_status);
meta("product.git_status_short", $git_short);

open my $ifh, "<", $manifest or die "$manifest: $!\n";
my $header = <$ifh>;
while (my $line = <$ifh>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    my ($name, $path, $region, $shape) = split /\t/, $line;
    next unless defined $name && length $name;
    next if %input_filter && !$input_filter{$name};
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
            my @cmd = test_view_cmd($mode, $threads, $path, $region);
            my $cmd_string = cmd_string(@cmd);
            for my $qsize (@qsizes) {
                for (1 .. $warmups) {
                    my ($status, undef, undef, $diag) = run_timed_loop($qsize, @cmd);
                    if ($status ne "ok") {
                        print STDERR "warmup failed: qsize=$qsize $cmd_string\n$diag\n";
                    }
                }
            }
            for my $repeat (1 .. $repeats) {
                my @rotated_qsizes = @qsizes;
                if (@rotated_qsizes) {
                    my $offset = ($repeat - 1) % @rotated_qsizes;
                    @rotated_qsizes = (@rotated_qsizes[$offset .. $#rotated_qsizes],
                                       @rotated_qsizes[0 .. $offset - 1]);
                }
                for my $qsize (@rotated_qsizes) {
                    my ($status, $seconds, $max_rss, $diag) = run_timed_loop($qsize, @cmd);
                    print $tfh join("\t", ts(), $mode, $threads, $qsize, $name,
                                    $shape, $region, $repeat, $loops,
                                    sprintf("%.6f", $seconds), $max_rss, $status,
                                    $cmd_string), "\n";
                    if ($status ne "ok") {
                        print STDERR "failed: qsize=$qsize $cmd_string\n$diag\n";
                    }
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
