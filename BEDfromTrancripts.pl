#!/usr/bin/perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;

my ($chrom_start, $thick_start, $RGB);

sub print_error{
print<<END_OF_LINE
cat perl BEDfromTranscripts -\"Options\"
     Options: 
           -T                Transcripts: the file of transcripts that you want to view
           -P                (Y or N):    prints transcripts which don\'t have coordinates in gtf file (use if repeats are removed)
           -L                Location of gtf file with the exon coordinates
END_OF_LINE
}

print_error() and die "Please Retry\n" if (!@ARGV); 

my %opts = @ARGV;

sub check_options{
    print_error and print "missing -T option, or file with transcripts\n" and die  unless (exists ($opts{-T}));
    print_error and print "missing -T option, or gtf file              \n"and die unless (exists ($opts{-L}));
}

check_options();

open my $transcripts,         $opts{-T} or die "Could not open $ARGV{-T}";
my $exon_coordinates =  $opts{-L};


while(<$transcripts>){
    chomp;
    my (@samples, @block_size, @block_start, @sample_group);
    my ($transcript_id, $XLOC, $gene, $cc, $rRNA,@genes) = split(/\t/);
    foreach my $g (@genes){
	next if ($g eq "-");
	my ($transcript) = split(/\|/,$g); # obtain 1st value of cufflinks column which should contain sample ID
	my ($sample, $sample_group);
	if ($transcript  =~m/q\d{1,2}:((CT|AD|D)\d{1,2})\.\d+/){ # grab sample  id information
	   $sample       = $1;
	   $sample_group = $2;
	}else{
	   $sample       = 'gene';
	   $sample_group = 'none';
	}

	push (@samples, $sample);
	push (@sample_group, $sample_group);
	$RGB = find_RGB(@sample_group);  # find RGB 
    }
    my $string = `grep $transcript_id  $exon_coordinates`;

    if ($string ne ''){
	my (@lines)   = split(/\n/,$string);
	my $new_gene  ='FALSE';
	my ($old_stop, $old_start, $thick_end, $chrom_stop, $chr, $cuff, $exon, $start, $stop, $strand, $r);
	$thick_start  = $chrom_start = $old_stop = $old_start = $thick_end = $chrom_stop = 0;
	my $gene_name      = join ('_',$transcript_id,@samples);

	foreach my $site (@lines){

	    ($chr, $cuff, $exon, $start, $stop, $r, $strand) = split(/\t/,$site);
	    find_start($start) if ($new_gene eq 'FALSE');

	     ###create exon block size: finds sizes of exons
	    my $block_size = ($stop - $start);        
	    my $block_start = ($start - $chrom_start); 
	    push ( @block_size, $block_size );
	    push ( @block_start, $block_start );
	    
	    ###make the old start(etc) the recent start so that if there are more exons, it will start from this recent exon
	    $old_start = $start;
	    $old_stop  = $stop;
	    $new_gene  = 'TRUE';
	    }

	(scalar(@block_size) == scalar (@block_start))? my $block_count = scalar (@block_size) : die "block arrays off"; # make exon blocks are equal!

	$thick_end = $chrom_stop = $old_stop;  # denote which areas should be filled in as an exon
	my $block_sizes = join(',',@block_size);
	my $block_start = join (',',@block_start);
	print join("\t", $chr, $chrom_start, $chrom_stop, $gene_name, "500", $strand, $thick_start,$thick_end, $RGB, $block_count, $block_sizes, $block_start)."\n";

    }
    else{
	print "$_\n" if ( (exists ($opts{-T})) and ($opts{-T}=~m/TRUE|T/) );
    }
}

sub find_RGB{
    my @samples = @_;
    my @uniq = uniq (@samples);
    my $sum = 0;
    foreach my $sample_type (@samples){
	my $n = $sample_type =~m/CT/i ? 4
	      : $sample_type =~m/AD/i ? 2
              :                         8;
	$sum += $n;
    }
    my $rgb = $sum==2  ? '255,0,0'    #red    if present in AD only
	    : $sum==4  ? '0,255,0'    #lime   if present in CT only 
	    : $sum==6  ? '0,0,255'    #blue   if present in AD + CT only
	    : $sum==8  ? '255,255,0'  #yellow if present in DLB only
	    : $sum==10 ? '128,0,128'  #purple if present in AD + DLB 
	    : $sum==12 ? '0,128,128'  #teal   if present in DLB + CT
	    :            '0,0,0';     #black  if present in all three groups
    return ($rgb);
}

sub find_start{
    my $start_position = shift;
    $chrom_start = $start_position;
    $thick_start = $start_position;
    return ($chrom_start, $thick_start);
}
    
