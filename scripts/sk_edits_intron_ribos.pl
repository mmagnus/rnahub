#!/usr/bin/perl -w
# sk_edits_intron_ribos.pl


use strict;
use Class::Struct;

use vars qw ($opt_v );  # required if strict used
use Getopt::Std;
getopts ('v');

# Print a helpful message if the user provides no input file.
if (@ARGV) {
        print "usage:  sk_edits_intron_ribos.pl [options]  \n\n";
        print "options:\n";
 	print "-v    :  be verbose\n";
 	exit;
}

my $maindir = "/Users/erivas/projects/SKennedy/2024_conserved_introns";
my $file_sk = "$maindir/data/sk_edits_intron_ribos.csv";
print "$file_sk\n";

my $file_sorek = "$maindir/data/sorek_suppl_info.csv";
print "$file_sk\n";

my $easel    = "~/src/Mysrc/R-scape/lib/hmmer/easel/miniapps/";
my $hmmbuild = "~/src/Mysrc/R-scape/lib/hmmer/src/hmmbuild";
my $nhmmer   = "~/src/Mysrc/R-scape/lib/hmmer/src/nhmmer";
my $hg38_fna = "$maindir/data/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna";
my $hg38_gff = "$maindir/data/ncbi_dataset/data/GCF_000001405.40/genomic.gff";

# SK 242, misses one entry: F10867_28
#
my $nseq_sk = 0;
my @entry_sk;
open (FILE, "$file_sk") || die;
while(<FILE>) {
    if (/^exon_id,/) {
	$entry_sk[$nseq_sk] = "";
    }
    elsif(/^(.+)\$/) {
	$entry_sk[$nseq_sk] .= $1;
	$nseq_sk ++;
	$entry_sk[$nseq_sk] = "";
    }
    else {
	$entry_sk[$nseq_sk] .= $_;
    }
    
}
close (FILE);
$nseq_sk ++;

print "nseq_sk = $nseq_sk\n";
if (0) {
    for (my $s = 0; $s < $nseq_sk; $s ++) {
	printf "%d\n$entry_sk[$s]\n\n", $s+1;
    }
}

# sorek 243
#-------------------------------------------------------------------------
my $nseq_sorek = 0;
my @entry_sorek;
open (FILE, "$file_sorek") || die;
while(<FILE>) {
    if (/Exon/) {
	$entry_sorek[$nseq_sorek] = "";
    }
    elsif(/^(.+)\n$/) {
	$entry_sorek[$nseq_sorek] .= $1;
	$nseq_sorek ++;
	$entry_sorek[$nseq_sorek] = "";
    }
    
}
close (FILE);

my $exon_filename               = "$maindir/data/sorek_suppl_info.exon";
my $intron_exon_intron_filename = "$maindir/data/sorek_suppl_info.intronU_exon_intronD";
my @exon_name;
parse_sorek($nseq_sorek, \@entry_sorek, \@exon_name, $exon_filename, $intron_exon_intron_filename);

#   create the 4 files with
#
#   Upstream_exon   U_intron     exon     D_intron     Downstream_exon
#   UUUUUUUUUUUUU---------------EEEEEE-----------------DDDDDDDDDDDDDDDDDD
#
#            xxxxxxxxxxxxxxxxxxxxxxxxx                              UE_UI_exon          
#                               xxxxxxxxxxxxxxxxxxxxxxxxxxxxx             exon_DI_DE
#                xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx               UI_exon_DI
#            xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx      UE_UI_exon_DI_DE
#
#            ----                 = flank_exon_chunk
#                 --------------- = intron chunk
my $flank_exon_chunk = 150;
my $intron_chunk     = 300;

my $file_exon             = "$maindir/data/sorek/sorek.exon.fa";
my $file_UE_UI_exon       = "$maindir/data/sorek/sorek.UE_UI_exon.fa";
my $file_exon_DI_DE       = "$maindir/data/sorek/sorek.exon_DI_DE.fa";
my $file_UI_exon_DI       = "$maindir/data/sorek/sorek.UI_exon_DI.fa";
my $file_UE_UI_exon_DI_DE = "$maindir/data/sorek/sorek.UE_UI_exon_DI_DE.fa";

open (my $F,    ">$file_exon")             || die;
open (my $FL,   ">$file_UE_UI_exon")       || die;
open (my $FR,   ">$file_exon_DI_DE")       || die;
open (my $FLR,  ">$file_UE_UI_exon_DI_DE") || die;
open (my $FLRI, ">$file_UI_exon_DI")       || die;

my $not_found = 0;
for (my $s = 0; $s < $nseq_sorek; $s ++) {

    # use sorek file with: 100nt_upstream+exon+100nt_downstream
    #                      to find the alt spliced exon in the annotation and extract the two adjacent exons and introns
    #
    my $filename = "$intron_exon_intron_filename.$exon_name[$s]";
    printf "\nNSEQ %d | $filename\n", $s+1;

    my $tblfile = run_hmmer($filename);

    my $this_file_exon             = "$maindir/data/sorek/sorek.exon/$exon_name[$s].exon.fa";
    my $this_file_UE_UI_exon       = "$maindir/data/sorek/sorek.UE_UI_exon/$exon_name[$s].UE_UI_exon.fa";
    my $this_file_exon_DI_DE       = "$maindir/data/sorek/sorek.exon_DI_DE/$exon_name[$s].exon_DI_DE.fa";
    my $this_file_UI_exon_DI       = "$maindir/data/sorek/sorek.UI_exon_DI/$exon_name[$s].UI_exon_DI.fa";
    my $this_file_UE_UI_exon_DI_DE = "$maindir/data/sorek/sorek.UE_UI_exon_DI_DE/$exon_name[$s].UE_UI_exon_DI_DE.fa";
    queryfiles_add($exon_name[$s], $tblfile, $hg38_gff, \$not_found, 
		   $F, $FL, $FR, $FLR, $FLRI,
		   $this_file_exon, $this_file_UE_UI_exon, $this_file_exon_DI_DE, $this_file_UI_exon_DI, $this_file_UE_UI_exon_DI_DE);
}
close($F);
close($FL);
close($FR);
close($FLR);
close($FLRI);
print "NOT FOUND $not_found\n";

sub parse_sorek {
    my ($nseq_sorek, $entry_sorek_ref, $exon_name_ref, $exon_filename, $intron_exon_intron_filename) = @_;

    my $exon_file               = "$exon_filename.fa";
    my $intron_exon_intron_file = "$intron_exon_intron_filename.fa";
    
    open (EF,    ">$exon_file") || die;
    open (IEIF, ">$intron_exon_intron_file") || die;

    print "nseq_sorek = $nseq_sorek\n";
    for (my $s = 0; $s < $nseq_sorek; $s ++) {
	#printf "\n%d\n$entry_sorek_ref->[$s]\n\n", $s+1;

	my @entry = split(/,/,$entry_sorek_ref->[$s]);
	my $ne = $#entry + 1;

	my $name     = $entry[0]; $exon_name_ref->[$s] = $name;
	my $U_intron = $entry[14];
	my $D_intron = $entry[16];
	my $exon     = $entry[12];
	my $desc     = "";
	for (my $c = 17; $c < $ne; $c ++) {
	    $desc  .= $entry[$c];
	}

	my $this_exon_file               = "$exon_filename.$name.fa";
	my $this_intron_exon_intron_file = "$intron_exon_intron_filename.$name.fa";
	open (TEF,   ">$this_exon_file") || die;
	open (TIEIF, ">$this_intron_exon_intron_file") || die;

	printf("\n$s-$ne> %s %s\n%s\n%s\n%s\n", $name, $desc, $U_intron, $exon, $D_intron);
	printf(EF    ">%s %s\n%s\n",     $name, $desc, $exon);
	printf(IEIF  ">%s %s\n%s%s%s\n", $name, $desc, $U_intron, $exon, $D_intron);
	printf(TEF   ">%s %s\n%s\n",     $name, $desc, $exon);
	printf(TIEIF ">%s %s\n%s%s%s\n", $name, $desc, $U_intron, $exon, $D_intron);

	close(TEF);
	close(TIEIF);

    }

    close(EF);
    close(IEIF);
}

sub queryfiles_add {
    my ($exon_name, $tblfile, $hg38_gff, $ret_not_found, 
	$F, $FL, $FR, $FLR, $FLRI, 
	$this_file_exon, $this_file_UE_UI_exon, $this_file_exon_DI_DE, $this_file_UI_exon_DI, $this_file_UE_UI_exon_DI_DE) = @_;
    
    my $exon_name_alt    = "";
    my $exon_strand_alt  = "";
    my $exon_from_alt    = -1;
    my $exon_to_alt      = -1;
    #upstream exon
    my $exon_name_Up     = "";
    my $exon_strand_Up   = "";
    my $exon_from_Up     = -1;
    my $exon_to_Up       = -1;
    #downstream exon
    my $exon_name_Down   = "";
    my $exon_strand_Down = "";
    my $exon_from_Down   = -1;
    my $exon_to_Down     = -1;
    
    my $not_found = find_updown_exons($tblfile, $hg38_gff, 
				      \$exon_name_alt,  \$exon_strand_alt,  \$exon_from_alt,  \$exon_to_alt,
				      \$exon_name_Up,   \$exon_strand_Up,   \$exon_from_Up,   \$exon_to_Up,
				      \$exon_name_Down, \$exon_strand_Down, \$exon_from_Down, \$exon_to_Down);
    $$ret_not_found += $not_found;
    
    # fetch the sequences
    # set the coordenates
    my $exon_len = $exon_to_alt - $exon_from_alt + 1;
    
    my $intronU_to   = $exon_from_alt  - 1;
    my $intronD_from = $exon_to_alt    + 1;
    my $intronU_from = ($exon_to_Up     > 0)? $exon_to_Up     + 1 : $intronU_to   - $intron_chunk;
    my $intronD_to   = ($exon_from_Down > 0)? $exon_from_Down - 1 : $intronD_from + $intron_chunk;
    my $intronU_len  = $intronU_to - $intronU_from + 1;
    my $intronD_len  = $intronD_to - $intronD_from + 1;

    my $frag_exonU_from  = ($exon_to_Up-$flank_exon_chunk > $exon_from_Up)? $exon_to_Up-$flank_exon_chunk : $exon_from_Up;
    my $frag_exonU_to    = $exon_to_Up;

    my $frag_exonD_from  = $exon_from_Down;
    my $frag_exonD_to    = ($exon_from_Down+$flank_exon_chunk < $exon_to_Down)? $exon_from_Down+$flank_exon_chunk : $exon_to_Down;

    my $frag_exonU_len   = ($frag_exonU_to > 0)? $frag_exonU_to - $frag_exonU_from + 1 : 0;
    my $frag_exonD_len   = ($frag_exonD_to > 0)? $frag_exonD_to - $frag_exonD_from + 1 : 0;

    print "not found? $not_found so far: $$ret_not_found\n";
    print "Up     exon: $exon_name_Up   $exon_strand_Down $frag_exonU_from $frag_exonU_to  len=$frag_exonU_len\n";
    print "Up   intron: $intronU_from   $intronU_to  len=$intronU_len\n";
    print "Alt    exon: $exon_name_alt  $exon_strand_alt  $exon_from_alt  $exon_to_alt  len=$exon_len\n";
    print "Down intron: $intronD_from   $intronD_to  len=$intronD_len\n";
    print "Down   exon: $exon_name_Down $exon_strand_Down $frag_exonD_from $frag_exonD_to  len=$frag_exonD_len\n";
    if ($intronU_len < 0 || $intronD_len < 0 || $frag_exonU_len < 0 || $frag_exonD_len < 0) {
	print "negative lengths!\n";
	die;
    }

    # fetch the sequences
    my $seqfile = "$tblfile.fa";
    my $cmd = "$easel/esl-sfetch $hg38_fna $exon_name_alt > $seqfile";
    system("$cmd\n");
    my $nsq = 0;
    my @seq;
    my @seq_name;
    parse_afafile($seqfile, \$nsq, \@seq, \@seq_name);
    if ($nsq != 1) { print "bad $seqfile\n"; die; }

    # introns are long, cap at $intron_chuck
    if ($intronU_len > $intron_chunk) {
	$intronU_from = $intronU_to - $intron_chunk + 1;
	$intronU_len  = $intronU_to - $intronU_from + 1;
	$frag_exonU_from = -1;
	$frag_exonU_to   = -1;
	$frag_exonU_len  = 0;
    }
    if ($intronD_len > $intron_chunk) {
	$intronD_to  = $intronD_from + $intron_chunk - 1;
	$intronD_len = $intronD_to - $intronD_from + 1;
	$frag_exonD_from = -1;
	$frag_exonD_to   = -1;
	$frag_exonD_len  = 0;
    }
    
    my $exon       = substr($seq[0], $exon_from_alt,   $exon_len);
    my $frag_exonU = substr($seq[0], $frag_exonU_from, $frag_exonU_len);
    my $frag_exonD = substr($seq[0], $frag_exonD_from, $frag_exonD_len);
    my $intronU    = substr($seq[0], $intronU_from,    $intronU_len);
    my $intronD    = substr($seq[0], $intronD_from,    $intronD_len);

    # set the names/sequences
    my $seq     =                    "$exon";
    my $seq_L   = "$frag_exonU$intronU$exon";
    my $seq_R   =                    "$exon$intronD$frag_exonD";
    my $seq_LRI =            "$intronU$exon$intronD";
    my $seq_LR  = "$frag_exonU$intronU$exon$intronD$frag_exonD";

    my $exon_coord    = $exon_from_alt."-".$exon_to_alt;
    my $intronU_coord = $intronU_from."-".$intronU_to;
    my $intronD_coord = $intronD_from."-".$intronD_to;
    my $frag_exonU_coord = ($frag_exonU_from > 0)? $frag_exonU_from."-".$frag_exonU_to : "0-0";
    my $frag_exonD_coord = ($frag_exonD_from > 0)? $frag_exonD_from."-".$frag_exonD_to : "0-0";
    
    my $seq_name     = $exon_name.".".$exon_name_alt.".".$exon_coord;
    my $seq_L_name   = $exon_name.".".$exon_name_alt.".".$frag_exonU_coord."/".$intronU_coord."/".$exon_coord;
    my $seq_R_name   = $exon_name.".".$exon_name_alt.".".$exon_coord."/".$intronD_coord."/".$frag_exonD_coord;
    my $seq_LRI_name = $exon_name.".".$exon_name_alt.".".$intronU_coord."/".$exon_coord."/".$intronD_coord;
    my $seq_LR_name  = $exon_name.".".$exon_name_alt.".".$frag_exonU_coord."/".$intronU_coord."/".$exon_coord."/".$intronD_coord."/".$frag_exonD_coord;
    print("exon:$seq_name\n");
    print("L:   $seq_L_name\n");
    print("R:   $seq_R_name\n");
    print("LRI: $seq_LRI_name\n");
    print("LR:  $seq_LR_name\n");

    printf $F    ">%s\n%s\n", $seq_name,     $seq;
    printf $FL   ">%s\n%s\n", $seq_L_name,   $seq_L;
    printf $FR   ">%s\n%s\n", $seq_R_name,   $seq_R;
    printf $FLRI ">%s\n%s\n", $seq_LRI_name, $seq_LRI;
    printf $FLR  ">%s\n%s\n", $seq_LR_name,  $seq_LR;
    
    open (my $THISF,    ">$this_file_exon")             || die;
    open (my $THISFL,   ">$this_file_UE_UI_exon")       || die;
    open (my $THISFR,   ">$this_file_exon_DI_DE")       || die;
    open (my $THISFLR,  ">$this_file_UE_UI_exon_DI_DE") || die;
    open (my $THISFLRI, ">$this_file_UI_exon_DI")       || die;
    
    printf $THISF    ">%s\n%s\n", $seq_name,     $seq;
    printf $THISFL   ">%s\n%s\n", $seq_L_name,   $seq_L;
    printf $THISFR   ">%s\n%s\n", $seq_R_name,   $seq_R;
    printf $THISFLRI ">%s\n%s\n", $seq_LRI_name, $seq_LRI;
    printf $THISFLR  ">%s\n%s\n", $seq_LR_name,  $seq_LR;

    close($THISF);
    close($THISFL);
    close($THISFR);
    close($THISFLRI);
    close($THISFLR);
    
    system("rm $seqfile\n");
}

sub find_updown_exons {
    my ($tblfile, $gfffile, 
	$ret_exon_name_alt,  $ret_exon_strand_alt,  $ret_exon_from_alt,  $ret_exon_to_alt,
	$ret_exon_name_Up,   $ret_exon_strand_Up,   $ret_exon_from_Up,   $ret_exon_to_Up,
	$ret_exon_name_Down, $ret_exon_strand_Down, $ret_exon_from_Down, $ret_exon_to_Down
	) = @_;

    my $exon_name   = "";
    my $exon_strand = "";
    my $exon_from   = -1;
    my $exon_to     = -1;
    parse_hmm_tbl($tblfile, \$exon_name, \$exon_from, \$exon_to, \$exon_strand);
    if ($exon_from >= $exon_to) { print "from = $exon_from shouls be < than to = $exon_to\ntblfile:$tblfile\n"; exit(); }
    print("^^hmm hit: $exon_name $exon_strand $exon_from $exon_to\n");

    my $exon_name_alt        = "";
    my $exon_strand_alt      = "";
    my $exon_from_alt        = -1;
    my $exon_to_alt          = -1;
    #upstream exon
    my $exon_name_Up     = "";
    my $exon_strand_Up   = "";
    my $exon_from_Up     = -1;
    my $exon_to_Up       = -1;
    #downstream exon
    my $exon_name_Down   = "";
    my $exon_strand_Down = "";
    my $exon_from_Down   = -1;
    my $exon_to_Down     = -1;
    
    my $gff_exon_name       = "";
    my $gff_exon_strand     = "";
    my $gff_exon_from       = -1;
    my $gff_exon_to         = -1;
    my $gff_exon_desc       = "";
    my $gff_exon_name_prv   = "";
    my $gff_exon_strand_prv = "";
    my $gff_exon_from_prv   = -1;
    my $gff_exon_to_prv     = -1;
    my $gff_exon_desc_prv   = -1;

    my $found     = 0;
    my $not_found = 0;
    my $n_exon = 0;
    open (FILE, "$gfffile") || die;
    while(<FILE>) {
	if (/\#/) {
	}
	elsif (/^($exon_name)\s+BestRefSeq\s+exon\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+/) {
	    $gff_exon_name   = $1;
	    $gff_exon_from   = $2;
	    $gff_exon_to     = $3;
	    $gff_exon_strand = $4;
	    $gff_exon_desc   = $5;

	    $n_exon ++;
	    
	    if ($exon_from <= $gff_exon_from && $exon_to >= $gff_exon_to) {
		$found ++;
		print "found exon! $found: $gff_exon_name $gff_exon_strand $gff_exon_from $gff_exon_to $gff_exon_desc\n";

		$exon_name_alt    = $gff_exon_name;
		$exon_from_alt    = $gff_exon_from;
		$exon_to_alt      = $gff_exon_to;
		$exon_strand_alt  = $gff_exon_strand;
		
		$exon_name_Up    = $gff_exon_name_prv;
		$exon_from_Up    = $gff_exon_from_prv;
		$exon_to_Up      = $gff_exon_to_prv;
		$exon_strand_Up  = $gff_exon_strand_prv;
		#print"exon # $n_exon\n\n";

	    }
	    elsif ($found > 0) {
		$exon_name_Down    = $gff_exon_name;
		$exon_from_Down    = $gff_exon_from;
		$exon_to_Down      = $gff_exon_to;
		$exon_strand_Down  = $gff_exon_strand;
		#print"next nexon # $n_exon\n";

		last;
	    }
	    
	    $gff_exon_name_prv   = $gff_exon_name;
	    $gff_exon_from_prv   = $gff_exon_from;
	    $gff_exon_to_prv     = $gff_exon_to;
	    $gff_exon_strand_prv = $gff_exon_strand;
	    $gff_exon_desc_prv   = $gff_exon_desc;
	}
    }
    close(FILE);

    # we haven't found this seq annotated as an exon, no exonU or exonD
    if ($found == 0) { 
	print "could not match to an exon hg38_gff\n"; 
	$not_found = 1;

	$exon_name_alt    = $exon_name;
	$exon_from_alt    = $exon_from + 100;
	$exon_to_alt      = $exon_to   - 100;
	$exon_strand_alt  = $exon_strand;
    }

    #print "Up    exon: $exon_name_Up   $exon_strand_Up   $exon_from_Up  $exon_to_Up\n";
    #print "Alt   exon: $exon_name_alt  $exon_strand_alt  $exon_from_alt  $exon_to_alt\n";
    #print "Down  exon: $exon_name_Down $exon_strand_Down $exon_from_Down $exon_to_Down\n";

    # reorder if necessary
    if ($exon_from_Up > $exon_to_Down) {
	my $foo_name   = $exon_name_Down;
	my $foo_from   = $exon_from_Down;
	my $foo_to     = $exon_to_Down;
	my $foo_strand = $exon_strand_Down;
	$exon_name_Down    = $exon_name_Up;
	$exon_from_Down    = $exon_from_Up;
	$exon_to_Down      = $exon_to_Up;
	$exon_strand_Down  = $exon_strand_Up;
	
	$exon_name_Up      = $foo_name;
	$exon_from_Up      = $foo_from;
	$exon_to_Up        = $foo_to;
	$exon_strand_Up    = $foo_strand;
    }

    $$ret_exon_name_alt    = $exon_name_alt;
    $$ret_exon_from_alt    = $exon_from_alt;
    $$ret_exon_to_alt      = $exon_to_alt;
    $$ret_exon_strand_alt  = $exon_strand_alt;

    $$ret_exon_name_Down   = $exon_name_Down;
    $$ret_exon_from_Down   = $exon_from_Down;
    $$ret_exon_to_Down     = $exon_to_Down;
    $$ret_exon_strand_Down = $exon_strand_Down;

    $$ret_exon_name_Up     = $exon_name_Up;
    $$ret_exon_from_Up     = $exon_from_Up;
    $$ret_exon_to_Up       = $exon_to_Up;
    $$ret_exon_strand_Up   = $exon_strand_Up;

    return $not_found;
}

sub run_hmmer {
    my ($filename) = @_;

    my $eval_thresh = 1e-3;
    
    my $fafile  = "$filename.fa";
    my $tblfile = "$filename.hmm.tbl";
    my $hmmfile = "$filename.hmm";

    if (-e $tblfile) {}
    else {
	my $cmd = "$hmmbuild $hmmfile $fafile";
	system("$cmd\n");
	
	$cmd = "$nhmmer --tblout $tblfile --incE $eval_thresh $hmmfile $hg38_fna";
	system("$cmd\n");
	
	system("rm $hmmfile\n");
    }
    
    return $tblfile;
}

sub parse_hmm_tbl {
    my ($tblfile, $ret_name, $ret_from, $ret_to, $ret_strand) = @_;

    my $hit = 0;
    my $name   = "";
    my $strand = "";
    my $from   = -1;
    my $to     = -1;
    open (FILE, "$tblfile") || die;
    while(<FILE>) {
	if (/\#/) {
	}
	elsif (/^(\S+)\s+\S\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+(\S)\s+/) {
	    $name   = $1;
	    $from   = $2;
	    $to     = $3;
	    $strand = $4;
	    $hit ++;
	    last;
	}
    }
    close(FILE);
    if ($hit == 0) { print "could not find hit in file $tblfile\n"; die; }

    if ($strand =~ "-") {
	my $ex = $from;
	$from  = $to;
	$to    = $ex
    }
    #print("^^hit $name $from $to\n");

    $$ret_name   = $name;
    $$ret_strand = $strand;
    $$ret_from   = $from;
    $$ret_to     = $to;
}

sub parse_afafile {
    my ($afafile, $ret_nsq, $asq_ref, $asqname_ref) = @_;

    my $n = 0;
    open(FILE, "$afafile") || die;
    while (<FILE>) {
	if (/>(\S+)\s*/) {
	    $asqname_ref->[$n] = $1;
	    $asq_ref->[$n]     = "";
	    $n ++;
	}
	elsif (/^(\S+)\s*$/) {
	    $asq_ref->[$n-1] .= $1;
	}
    }
    close(FILE);

    my $maxlen = 0;
    for (my $x = 0; $x < $n; $x ++) {
	if (length($asqname_ref->[$x]) > $maxlen) { $maxlen = length($asqname_ref->[$x]); }
    }
    
    for (my $x = 0; $x < $n; $x ++) {
	while (length($asqname_ref->[$x]) < $maxlen) { $asqname_ref->[$x] .= " "; }
    }
    
    $$ret_nsq = $n;
    return length($asq_ref->[0]);
}
