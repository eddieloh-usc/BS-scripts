#!/usr/bin/perl

use strict;

my $configfile = $ARGV[0] or die "Require path to configuration file\n"; 

my $starttime=localtime();
print STDERR "\n>>> Analysis started on $starttime\n\n";

##### Parse Config File #####
print STDERR ">>> Parsing configuration file...\n";
my %config=();
my $countconfigline=0;
open (CONFIG,"$configfile") or die "ERROR! Cannot open configuration file $configfile!!!\n";
while (my $line=<CONFIG>) {
   chomp $line;
   $line=~s/^\s+//gm;
   $line=~s/\s+$//gm;
   $line=~s/#.*//gm;
   my $line2=$line;
   next if !($line =~ /\S/gm);
   my @array=split(/\s+/,$line2);
   my $arraylen=scalar(@array);
   next if ($array[0] eq "Parameter"); #skip header line
   die "ERROR! $arraylen columns parsed where only 2 columns expected!!!\n" if ($arraylen != 2);
   $countconfigline++;
   $config{$array[0]}=$array[1];
   print STDERR "$countconfigline   $array[0]    \t$array[1]\n";
}
my $numparam=scalar(keys(%config));
print STDERR ">>> Configuration file parsed... $numparam parameter(s) stored...\n";
#sleep(5);

##### Check All Parameters Present and Load #####
die "Fastq_Dir parameter not provided\n\n" if (!exists $config{"Fastq_Dir"});
die "Output_Dir parameter not provided\n\n" if (!exists $config{"Output_Dir"});
my $fastqdir = $config{"Fastq_Dir"};
my $outputdir = $config{"Output_Dir"};

die "Insert_File parameter not provided\n\n" if (!exists $config{"Insert_File"});
die "Left_Adapter parameter not provided\n\n" if (!exists $config{"Left_Adapter"});
die "Right_Adapter parameter not provided\n\n" if (!exists $config{"Right_Adapter"});
die "Trim_Length parameter not provided\n\n" if (!exists $config{"Trim_Length"});
my $insertfile = $config{"Insert_File"};
my $leftadapter = $config{"Left_Adapter"};
my $rightadapter = $config{"Right_Adapter"};
my $trimlength = $config{"Trim_Length"};

die "Run_Trim parameter not provided\n\n" if (!exists $config{"Run_Trim"});
die "Remove_dup parameter not provided\n\n" if (!exists $config{"Remove_dup"});
my $runtrim = $config{"Run_Trim"};
my $removedup = $config{"Remove_dup"};

die "Quality_cutoff parameter not provided\n\n" if (!exists $config{"Quality_cutoff"});
die "R1match_cutoff parameter not provided\n\n" if (!exists $config{"R1match_cutoff"});
die "R2match_cutoff parameter not provided\n\n" if (!exists $config{"R2match_cutoff"});
die "R1others_cutoff parameter not provided\n\n" if (!exists $config{"R1others_cutoff"});
die "R2others_cutoff parameter not provided\n\n" if (!exists $config{"R2others_cutoff"});
die "R1ct_cutoff_A parameter not provided\n\n" if (!exists $config{"R1ct_cutoff_A"});
die "R2ct_cutoff_A parameter not provided\n\n" if (!exists $config{"R2ct_cutoff_A"});
die "R1ga_cutoff_A parameter not provided\n\n" if (!exists $config{"R1ga_cutoff_A"});
die "R2ga_cutoff_A parameter not provided\n\n" if (!exists $config{"R2ga_cutoff_A"});
die "R1ga_cutoff_B parameter not provided\n\n" if (!exists $config{"R1ga_cutoff_B"});
die "R2ga_cutoff_B parameter not provided\n\n" if (!exists $config{"R2ga_cutoff_B"});
die "R1ct_cutoff_B parameter not provided\n\n" if (!exists $config{"R1ct_cutoff_B"});
die "R2ct_cutoff_B parameter not provided\n\n" if (!exists $config{"R2ct_cutoff_B"});
my $qualcut = $config{"Quality_cutoff"};
my $R1matchcut = $config{"R1match_cutoff"}; 
my $R2matchcut = $config{"R2match_cutoff"};
my $R1otherscut = $config{"R1others_cutoff"};
my $R2otherscut = $config{"R2others_cutoff"};
my $R1ctcutA = $config{"R1ct_cutoff_A"};
my $R2ctcutA = $config{"R2ct_cutoff_A"};
my $R1gacutA = $config{"R1ga_cutoff_A"};
my $R2gacutA = $config{"R2ga_cutoff_A"};
my $R1gacutB = $config{"R1ga_cutoff_B"};
my $R2gacutB = $config{"R2ga_cutoff_B"};
my $R1ctcutB = $config{"R1ct_cutoff_B"};
my $R2ctcutB = $config{"R2ct_cutoff_B"};

#die "Script04 parameter not provided\n\n" if (!exists $config{"R_PlotScript"});
#die "Plot_Highlight parameter not provided\n\n" if (!exists $config{"Plot_Highlight"});
#die "Plot_ymin parameter not provided\n\n" if (!exists $config{"Plot_ymin"});
#die "Plot_ymax parameter not provided\n\n" if (!exists $config{"Plot_ymax"});
#die "Plot_xmin parameter not provided\n\n" if (!exists $config{"Plot_xmin"});
#die "Plot_xmax parameter not provided\n\n" if (!exists $config{"Plot_xmax"});
#my $Rplotscript = $config{"R_PlotScript"};
#my $hotstring = $config{"Plot_Highlight"};
#my $plotymin = $config{"Plot_ymin"};
#my $plotymax = $config{"Plot_ymax"};
#my $plotxmin = $config{"Plot_xmin"};
#my $plotxmax = $config{"Plot_xmax"};

print STDERR ">>> All required  parameter(s) present and loaded...\n\n";

##### Check Executables Exist #####
print STDERR ">>> Checking progam executables and scripts...\n";
use Env '@PATH';
my $prog2="cutadapt";
my $prog2exist = grep -x "$_/$prog2", @PATH;
die "cutadapt program not found in PATH, please check!!!\n" if ($prog2exist != 1);
print STDERR ">>> Program check completed...\n\n";

##### Make Output directory #####
print STDERR ">>> Checking Output Directory...\n";
if (-e $outputdir) {
   #do nothing
}else {
   system("mkdir $outputdir");
   print STDERR "Output directory $outputdir does not exist and is created...\n";
}
print STDERR ">>> Output directory check completed...\n\n";

##### Make Log directory #####
print STDERR ">>> Checking Log Directory...\n";
my $logdir="$outputdir/RunLogs";
if (-e $logdir) {
   #do nothing
}else {
   system("mkdir $logdir");
   print STDERR "Log directory $logdir does not exist and is created...\n";
}
print STDERR ">>> Log directory check completed...\n\n";

##### Make Trim directory #####
print STDERR ">>> Checking Trimmed Sequences Output Directory...\n";
my $trimmeddir="$outputdir/TrimmedFastq";
if (-e $trimmeddir) {
   #do nothing
}else {
   system("mkdir $trimmeddir");
   print STDERR "Trimmed directory $trimmeddir does not exist and is created...\n";
}
print STDERR ">>> Trimmed directory check completed...\n\n";


##### Define Functions #####
### reverse complement ###
sub revcomp {    
  my ($forward) = @_;
  my @forward=split(//,$forward);
  my $rc="";
  while (my $base=pop(@forward)) {
     if ($base eq "A") {
        $rc=$rc."T";
     }elsif ($base eq "T") {
        $rc=$rc."A";
     }elsif ($base eq "G") {
        $rc=$rc."C";
     }elsif ($base eq "C") {
        $rc=$rc."G";
     }elsif ($base eq "N") {
        $rc=$rc."N";
     }else {
        die "RevComp error! $forward...$base...\n";
     }
  }
  return($rc);
}
### reverse quality score ###
sub revqual {
  my ($forward) = @_;
  $forward =~ s/0/1/gm;  #some regex issue with 0's, ok to switch to 1's
  my @forward=split(//,$forward);
  my $forwardlen=scalar(@forward);
  my $rq="";
  while (my $basequal1=pop(@forward)) {
     $rq=$rq."$basequal1";
  }
  return($rq);
}
### bin quality scores into scale from 1-8 ###
sub binqual {
  my ($forward) = @_;
  $forward =~ s/0/1/gm;  #some regex issue with 0's, ok to switch to 1's
  my @forward=split(//,$forward);
  my $forwardlen=scalar(@forward);
  my $rq="";
  while (my $basequal1=shift(@forward)) {
     my $basequal2=$basequal1; my $basequal3=$basequal1;
     my $basequal4=$basequal1; my $basequal5=$basequal1;
     my $basequal6=$basequal1; my $basequal7=$basequal1;
     my $basequal8=$basequal1;
     if ($basequal1 =~ /[\!\"\#\$\%]/gm) {
        $rq=$rq."1";
     }elsif ($basequal2 =~ /[\&\'\(\)\*]/gm) {
        $rq=$rq."2";
     }elsif ($basequal3 =~ /[\+\,\-\.\/]/gm) {
        $rq=$rq."3";
     }elsif ($basequal4 =~ /[01234]/gm) {
        $rq=$rq."4";
     }elsif ($basequal5 =~ /[56789]/gm) {
        $rq=$rq."5";
     }elsif ($basequal6 =~ /[\:\;\<\=\>]/gm) {
        $rq=$rq."6";
     }elsif ($basequal7 =~ /[\?\@ABC]/gm) {
        $rq=$rq."7";
     }elsif ($basequal8 =~ /[DEFGHI]/gm) {
        $rq=$rq."8";
     }else {
        die "RevQual error! $forward...$basequal1...\n";
     }
  }
  return($rq);
}
### fillhash ###
sub fillhash {
   my ($hash,$string1,$string2) = @_;
   my %hash=%$hash;
   my $string=$string1.$string2;

   while ($string=~/([^\:]+)\:\:/gm) {
      my $temp=$1;
      next if !($temp=~/\,/);
      $temp=~s/pos//gm;
      $temp=~s/LOWQUAL/LQ/gm;
      #$temp=~s/^\S//gm;
      #print STDERR "  $temp\n";
      my @temp=split(/\,/,$temp);
      my $mut=$temp[0];
      for (my $i=1;$i<=scalar(@temp);$i++) {
         if (!exists $hash{$temp[$i]}) {
            $hash{$temp[$i]}=$mut;
         }elsif ($hash{$temp[$i]} eq $mut) {
            #do nothing 
         }elsif ($hash{$temp[$i]} ne $mut) {
            $hash{$temp[$i]}="$hash{$temp[$i]}\:\:$mut";
         }else {
            die "ERROR3!\n";
	 }
      }
   }
   foreach my $key (sort {$a <=> $b} (keys %hash)) {
      #print STDERR "\t$key\_$hash{$key}";
      if ($key eq "") {
	 #do nothing
      }elsif ($hash{$key}=~/^[ATGC][ATGC]$/gm) {
	 #do nothing
      }elsif ($hash{$key}=~/^LQ$/gm) {
	 #do nothing
      }elsif ($hash{$key}=~/LQ\:\:/gm) {
         $hash{$key}=~s/LQ\:\://gm;
      }elsif ($hash{$key}=~/\:\:LQ/gm) {
         $hash{$key}=~s/\:\:LQ//gm;
      }elsif ($hash{$key}=~/\:\:/gm) {
	 $hash{$key}="MM";
	 #print STDERR "  pause $hash{$key}  ";
	 #sleep 5;
      }else {
	 die "  ERROR4!\t$hash{$key}\n\n";
      }
      #print STDERR "\t->\t$key\_$hash{$key}\n";
   }
   return %hash;
}


##### Processing Begins #####
##### Trimming Fastq Files #####
print STDERR ">>> Fastq File Trimming...\n";
if ($runtrim eq "Yes") {
   my $countfastqfile=0;
   foreach my $file1 (glob "$fastqdir/*_R1_*fastq*") {
      $countfastqfile=$countfastqfile+2;
      my $file1out1=$file1;
      $file1out1=~s/$fastqdir/$trimmeddir/gm;
      $file1out1=~s/\.fastq/\.pretrimmed\.fastq/gm;
      my $file1out2=$file1out1;
      $file1out2=~s/pretrimmed/trimmed/gm;
      my $logfile=$file1out2;
      $logfile=~s/$trimmeddir/$logdir/gm;
      $logfile=~s/_R1_/_R1n2_/gm;
      $logfile=~s/trimmed.fastq\S*/01_trim.log/gm;

      my $file2=$file1;
      $file2=~s/\_R1\_/\_R2\_/gm;
      die "ERROR!!! Cannot find R2 file for $file1...\n\n" if !(-e $file2);
      my $file2out1=$file2;
      $file2out1=~s/$fastqdir/$trimmeddir/gm;
      $file2out1=~s/\.fastq/\.pretrimmed\.fastq/gm;
      my $file2out2=$file2out1;
      $file2out2=~s/pretrimmed/trimmed/gm;

      print STDERR "Trimming $file1\n";
      print STDERR "Trimming $file2\n";

      my $cmd1="cutadapt -g $leftadapter -G $rightadapter --pair-filter=any --discard-untrimmed -o $file1out1 -p $file2out1 $file1 $file2  > $logfile";
      my $cmd2="cutadapt -l $trimlength -m $trimlength --pair-filter=any -o $file1out2 -p $file2out2 $file1out1 $file2out1 >> $logfile";
      my $cmd3="rm $file1out1 $file2out1";
      system($cmd1);
      system($cmd2);
      system($cmd3);
      system("echo \" \" >> $logfile");

      if ($file1out2=~/gz$/gm) { system("gunzip $file1out2"); }
      if ($file2out2=~/gz$/gm) { system("gunzip $file2out2"); }
   }
}elsif ($runtrim eq "No") {
   print STDERR "Trimming skipped as set by Run_Trim in config file...\n";
}else {
   die "ERROR! Run_Trim needs to be set to \"Yes\" or \"No\" only. No other settings accepted...\n";
}  
print STDERR ">>> Trimming fastq files completed...\n\n";


##### Sequence Analysis #####
print STDERR ">>> Running Sequence Analysis...\n";
foreach my $inputfile1 (glob "$trimmeddir/*_R1_*fastq") {
   my $inputfile2 = $inputfile1;
   $inputfile2 =~ s/R1/R2/gm;

   my $sample=$inputfile1;
   $sample=~s/$trimmeddir/$outputdir/gm;
   $sample=~s/_R1_/_R1n2_/gm;
   $sample=~s/\.trimmed.fastq//gm;

   my $outputfile1=$sample.".01_allreads_mismatch_stats.txt";
   my $outputfile2=$sample.".02_aln_by_category.txt";
   my $logfile=$sample.".02_analysis.log";
   $logfile=~s/$outputdir/$logdir/gm;

   unless (open(LOGFILE, ">$logfile")){
      print STDERR "Cannot open file \"$logfile\" to write to!!\n\n";
      exit;
   }
   print STDERR "\nProcessing $inputfile1\n";
   print LOGFILE "\nProcessing $inputfile1\n";

   unless (open(OUTFILE, ">$outputfile1")){
      print STDERR "Cannot open file \"$outputfile1\" to write to!!\n\n";
      exit;
   }
   print OUTFILE "number\treadID\torientation\tnumMatch\tnumMismatch\tnumgoodCT\tnumgoodGA\tnumOthersPlusLowqual\treadSequence\tgoodATpositions\tgoodAGpositions\tgoodACpositions\tgoodTApositions\tgoodTGpositions\tgoodTCpositions\tgoodGApositions\tgoodGTpositions\tgoodGCpositions\tgoodCApositions\tgoodCTpositions\tgoodCGpositions\tLowQualpositions\tfilterCategory\n";
   my $PPout=$outputfile2;
   my $PFout=$outputfile2;
   my $FPout=$outputfile2;
   my $FFout=$outputfile2;
   $PPout=~s/\.txt$/\_Cond1Pass_Cond2Pass\.txt/gm;
   $PFout=~s/\.txt$/\_Cond1Pass_Cond2Fail\.txt/gm;
   $FPout=~s/\.txt$/\_Cond1Fail_Cond2Pass\.txt/gm;
   $FFout=~s/\.txt$/\_Cond1Fail_Cond2Fail\.txt/gm;
   die "Cannot open $PPout" unless (open(PP, ">$PPout"));
   die "Cannot open $PFout" unless (open(PF, ">$PFout"));
   die "Cannot open $FPout" unless (open(FP, ">$FPout"));
   die "Cannot open $FFout" unless (open(FF, ">$FFout"));

   ### Process reference sequence ###
   print STDERR "Part1 - Read reference sequence\n";
   print LOGFILE "Part1 - Read reference sequence\n";
   my $refid;
   my $refseq;
   open (INPUT,"$insertfile") or die "Cannot open insert reference file!!! \n";
   while (my $line=<INPUT>) {
      chomp $line;
      my $line2=$line;
      if ($line2 =~ />(\S+)/gm) {
         $refid=$1;
      }else {
         $line2 =~ s/\s//gm;
         $refseq=$refseq.$line2;
      }
   }
   $refseq=uc($refseq);
   my $reflen=length($refseq);
   my @refseq=split(//,$refseq);
   my $reflen2=scalar(@refseq);
   print STDERR ">$refid\t$reflen\t$reflen2\n$refseq\n";
   print LOGFILE ">$refid\t$reflen\t$reflen2\n$refseq\n";

   my $paddash="";
   my $padN="";
   for (my $i=1;$i<=($reflen-$trimlength);$i++) {
      $paddash=$paddash."-";
      $padN=$padN."N";
   }

   print STDERR "Part2 - Process fastq file\nParsing...";
   print LOGFILE "Part2 - Process fastq file\nParsing...";
   ### Initialize variables ###
   my $countread=0; my $countline=0; my $countpairs=0;
   my $readidR1=""; my $readidR2="";
   my $readseqR1=""; my $readseqR2="";
   my $nummatchR1=0; my $nummatchR2=0;
   my $nummismatchR1=0; my $nummismatchR2=0;
   my $numatR1=0; my $numatR2=0;
   my $numagR1=0; my $numagR2=0;
   my $numacR1=0; my $numacR2=0;
   my $numtaR1=0; my $numtaR2=0;
   my $numtgR1=0; my $numtgR2=0;
   my $numtcR1=0; my $numtcR2=0;
   my $numgaR1=0; my $numgaR2=0;
   my $numgtR1=0; my $numgtR2=0;
   my $numgcR1=0; my $numgcR2=0;
   my $numcaR1=0; my $numcaR2=0;
   my $numctR1=0; my $numctR2=0;
   my $numcgR1=0; my $numcgR2=0;
   my $numothersR1=0; my $numothersR2=0;
   my $numlqR1=0; my $numlqR2=0;
   my $atpositionsR1="ATpos"; my $atpositionsR2="ATpos";
   my $agpositionsR1="AGpos"; my $agpositionsR2="AGpos";
   my $acpositionsR1="ACpos"; my $acpositionsR2="ACpos";
   my $tapositionsR1="TApos"; my $tapositionsR2="TApos";
   my $tgpositionsR1="TGpos"; my $tgpositionsR2="TGpos";
   my $tcpositionsR1="TCpos"; my $tcpositionsR2="TCpos";
   my $gapositionsR1="GApos"; my $gapositionsR2="GApos";
   my $gtpositionsR1="GTpos"; my $gtpositionsR2="GTpos";
   my $gcpositionsR1="GCpos"; my $gcpositionsR2="GCpos";
   my $capositionsR1="CApos"; my $capositionsR2="CApos";
   my $ctpositionsR1="CTpos"; my $ctpositionsR2="CTpos";
   my $cgpositionsR1="CGpos"; my $cgpositionsR2="CGpos";
   my $otpositionsR1="OTHERpos"; my $otpositionsR2="OTHERpos";
   my $lqpositionsR1="LOWQUALpos"; my $lqpositionsR2="LOWQUALpos";
   my $alignseqR1=""; my $alignseqR2="";
   my $countPP=0;
   my $countPF=0;
   my $countFP=0;
   my $countFF=0;

   open (INPUT1,"$inputfile1") or die "Cannot open R1 input fastq file!!! \n";
   open (INPUT2,"$inputfile2") or die "Cannot open R2 input fastq file!!! \n";
   while ((my $lineR1=<INPUT1>)&&(my $lineR2=<INPUT2>)) {
      chomp $lineR1;
      my $line2R1=$lineR1;
      my $line3R1=$lineR1;
      chomp $lineR2;
      my $line2R2=$lineR2;
      my $line3R2=$lineR2;

      if ($line2R2 =~ /\@(M\S+)/gm) {
         $readidR2=$1;
         $line2R1 =~ /\@(M\S+)/gm or die "mismatch line...\n>>$line2R1\n>>$line2R2\n";
         $readidR1=$1;
         die "mismatch id...$readidR1...$readidR2...\n" if ($readidR1 ne $readidR2);

         $countread++;
         $countline=1;
      }else {
         $countline++;
      }

      if ($countline==2) {
         $countpairs++;
         $lineR1=~s/\s//gm;
         $lineR2=~s/\s//gm;

         $readseqR1=uc($lineR1);
         $readseqR1=$readseqR1.$padN;
         $readseqR2=revcomp(uc($lineR2));
         $readseqR2=$padN.$readseqR2;

	 if ($readseqR1=~/[^ATGCN]/gm) {
            die "unknown base in R1: $readseqR1\n";
         }elsif ($readseqR2=~/[^ATGCN]/gm) {
            die "unknown base in R2: $readseqR2\n";
         }
      }elsif ($countline==4) {
         $line3R1=~s/\s//gm;
         $line3R2=~s/\s//gm;

         my $qualR1=binqual($line3R1).$padN;
         my $qualR2=$padN.revqual(binqual($line3R2));

         # further process reads and quals
         my @readseqR1=split(//,$readseqR1);
         my @readseqR2=split(//,$readseqR2);
         my $readlenR1=scalar(@readseqR1);
         my $readlenR2=scalar(@readseqR2);
         die "ERROR! readlenR1 $readlenR1...\n\n" if ($readlenR1!=$reflen);
         die "ERROR! readlenR2 $readlenR2...\n\n" if ($readlenR2!=$reflen);

         my @qualR1=split(//,$qualR1);
         my @qualR2=split(//,$qualR2);
         my $quallenR1=scalar(@qualR1);
         my $quallenR2=scalar(@qualR2);
         die "ERROR! quallenR1 $quallenR1...$readidR1...$qualR2...\n\n" if ($quallenR1!=$reflen);
         die "ERROR! quallenR2 $quallenR2...$readidR2...$qualR2...\n\n" if ($quallenR2!=$reflen);

         #print STDERR "QUALCUT\t$qualcut\n";

	 #R1
         for (my $i=0;$i<$trimlength;$i++) {
            if ($refseq[$i] eq $readseqR1[$i]) {
               $nummatchR1++;
               $alignseqR1=$alignseqR1."\.";
            }else {
               $nummismatchR1++;
               if (($refseq[$i] eq "A")&&($readseqR1[$i] eq "T")&&($qualR1[$i] >= $qualcut)) {
                  $numatR1++;
                  $atpositionsR1=$atpositionsR1.",".($i+1);
                  $numothersR1++;
                  $otpositionsR1=$otpositionsR1.",".($i+1);
	       }elsif (($refseq[$i] eq "A")&&($readseqR1[$i] eq "G")&&($qualR1[$i] >= $qualcut)) {
                  $numagR1++;
                  $agpositionsR1=$agpositionsR1.",".($i+1);
                  $numothersR1++;
                  $otpositionsR1=$otpositionsR1.",".($i+1);
	       }elsif (($refseq[$i] eq "A")&&($readseqR1[$i] eq "C")&&($qualR1[$i] >= $qualcut)) {
                  $numacR1++;
                  $acpositionsR1=$acpositionsR1.",".($i+1);
                  $numothersR1++;
                  $otpositionsR1=$otpositionsR1.",".($i+1);
               }elsif (($refseq[$i] eq "T")&&($readseqR1[$i] eq "A")&&($qualR1[$i] >= $qualcut)) {
                  $numtaR1++;
                  $tapositionsR1=$tapositionsR1.",".($i+1);
                  $numothersR1++;
                  $otpositionsR1=$otpositionsR1.",".($i+1);
               }elsif (($refseq[$i] eq "T")&&($readseqR1[$i] eq "G")&&($qualR1[$i] >= $qualcut)) {
                  $numtgR1++;
                  $tgpositionsR1=$tgpositionsR1.",".($i+1);
                  $numothersR1++;
                  $otpositionsR1=$otpositionsR1.",".($i+1);
               }elsif (($refseq[$i] eq "T")&&($readseqR1[$i] eq "C")&&($qualR1[$i] >= $qualcut)) {
                  $numtcR1++;
                  $tcpositionsR1=$tcpositionsR1.",".($i+1);
                  $numothersR1++;
                  $otpositionsR1=$otpositionsR1.",".($i+1);
               }elsif (($refseq[$i] eq "G")&&($readseqR1[$i] eq "A")&&($qualR1[$i] >= $qualcut)) {
                  $numgaR1++;
                  $gapositionsR1=$gapositionsR1.",".($i+1);
               }elsif (($refseq[$i] eq "G")&&($readseqR1[$i] eq "T")&&($qualR1[$i] >= $qualcut)) {
                  $numgtR1++;
                  $gtpositionsR1=$gtpositionsR1.",".($i+1);
                  $numothersR1++;
                  $otpositionsR1=$otpositionsR1.",".($i+1);
               }elsif (($refseq[$i] eq "G")&&($readseqR1[$i] eq "C")&&($qualR1[$i] >= $qualcut)) {
                  $numgcR1++;
                  $gcpositionsR1=$gcpositionsR1.",".($i+1);
                  $numothersR1++;
                  $otpositionsR1=$otpositionsR1.",".($i+1);
               }elsif (($refseq[$i] eq "C")&&($readseqR1[$i] eq "A")&&($qualR1[$i] >= $qualcut)) {
                  $numcaR1++;
                  $capositionsR1=$capositionsR1.",".($i+1);
                  $numothersR1++;
                  $otpositionsR1=$otpositionsR1.",".($i+1);
               }elsif (($refseq[$i] eq "C")&&($readseqR1[$i] eq "T")&&($qualR1[$i] >= $qualcut)) {
                  $numctR1++;
                  $ctpositionsR1=$ctpositionsR1.",".($i+1);
               }elsif (($refseq[$i] eq "C")&&($readseqR1[$i] eq "G")&&($qualR1[$i] >= $qualcut)) {
                  $numcgR1++;
                  $cgpositionsR1=$cgpositionsR1.",".($i+1);
                  $numothersR1++;
                  $otpositionsR1=$otpositionsR1.",".($i+1);
               }elsif (($refseq[$i] ne $readseqR1[$i])&&($qualR1[$i] < $qualcut)) {
                  $numlqR1++;
                  $lqpositionsR1=$lqpositionsR1.",".($i+1);
               }else {
		  #$numothersR1++;
		  #$otpositionsR1=$otpositionsR1.",".($i+1);
		  die "unknown R1 base comparison: $refseq[$i]\_$readseqR1[$i]\_$qualR1[$i]\n";
               }
               $alignseqR1=$alignseqR1."$readseqR1[$i]";
            }
         }
         $alignseqR1=$alignseqR1.$padN;


	 #R2
         for (my $i=($reflen-$trimlength);$i<$reflen;$i++) {
            if ($refseq[$i] eq $readseqR2[$i]) {
               $nummatchR2++;
               $alignseqR2=$alignseqR2."\.";
            }else {
               $nummismatchR2++;
               if (($refseq[$i] eq "A")&&($readseqR2[$i] eq "T")&&($qualR2[$i] >= $qualcut)) {
                  $numatR2++;
                  $atpositionsR2=$atpositionsR2.",".($i+1);
                  $numothersR2++;
                  $otpositionsR2=$otpositionsR2.",".($i+1);
	       }elsif (($refseq[$i] eq "A")&&($readseqR2[$i] eq "G")&&($qualR2[$i] >= $qualcut)) {
                  $numagR2++;
                  $agpositionsR2=$agpositionsR2.",".($i+1);
                  $numothersR2++;
                  $otpositionsR2=$otpositionsR2.",".($i+1);
	       }elsif (($refseq[$i] eq "A")&&($readseqR2[$i] eq "C")&&($qualR2[$i] >= $qualcut)) {
                  $numacR2++;
                  $acpositionsR2=$acpositionsR2.",".($i+1);
                  $numothersR2++;
                  $otpositionsR2=$otpositionsR2.",".($i+1);
               }elsif (($refseq[$i] eq "T")&&($readseqR2[$i] eq "A")&&($qualR2[$i] >= $qualcut)) {
                  $numtaR2++;
                  $tapositionsR2=$tapositionsR2.",".($i+1);
                  $numothersR2++;
                  $otpositionsR2=$otpositionsR2.",".($i+1);
               }elsif (($refseq[$i] eq "T")&&($readseqR2[$i] eq "G")&&($qualR2[$i] >= $qualcut)) {
                  $numtgR2++;
                  $tgpositionsR2=$tgpositionsR2.",".($i+1);
                  $numothersR2++;
                  $otpositionsR2=$otpositionsR2.",".($i+1);
               }elsif (($refseq[$i] eq "T")&&($readseqR2[$i] eq "C")&&($qualR2[$i] >= $qualcut)) {
                  $numtcR2++;
                  $tcpositionsR2=$tcpositionsR2.",".($i+1);
                  $numothersR2++;
                  $otpositionsR2=$otpositionsR2.",".($i+1);
               }elsif (($refseq[$i] eq "G")&&($readseqR2[$i] eq "A")&&($qualR2[$i] >= $qualcut)) {
                  $numgaR2++;
                  $gapositionsR2=$gapositionsR2.",".($i+1);
               }elsif (($refseq[$i] eq "G")&&($readseqR2[$i] eq "T")&&($qualR2[$i] >= $qualcut)) {
                  $numgtR2++;
                  $gtpositionsR2=$gtpositionsR2.",".($i+1);
                  $numothersR2++;
                  $otpositionsR2=$otpositionsR2.",".($i+1);
               }elsif (($refseq[$i] eq "G")&&($readseqR2[$i] eq "C")&&($qualR2[$i] >= $qualcut)) {
                  $numgcR2++;
                  $gcpositionsR2=$gcpositionsR2.",".($i+1);
                  $numothersR2++;
                  $otpositionsR2=$otpositionsR2.",".($i+1);
               }elsif (($refseq[$i] eq "C")&&($readseqR2[$i] eq "A")&&($qualR2[$i] >= $qualcut)) {
                  $numcaR2++;
                  $capositionsR2=$capositionsR2.",".($i+1);
                  $numothersR2++;
                  $otpositionsR2=$otpositionsR2.",".($i+1);
               }elsif (($refseq[$i] eq "C")&&($readseqR2[$i] eq "T")&&($qualR2[$i] >= $qualcut)) {
                  $numctR2++;
                  $ctpositionsR2=$ctpositionsR2.",".($i+1);
               }elsif (($refseq[$i] eq "C")&&($readseqR2[$i] eq "G")&&($qualR2[$i] >= $qualcut)) {
                  $numcgR2++;
                  $cgpositionsR2=$cgpositionsR2.",".($i+1);
                  $numothersR2++;
                  $otpositionsR2=$otpositionsR2.",".($i+1);
               }elsif (($refseq[$i] ne $readseqR2[$i])&&($qualR2[$i] < $qualcut)) {
                  $numlqR2++;
                  $lqpositionsR2=$lqpositionsR2.",".($i+1);
               }else {
		  #$numothersR2++;
		  #$otpositionsR2=$otpositionsR2.",".($i+1);
		  die "unknown R2 base comparison: $refseq[$i]\_$readseqR2[$i]\_$qualR2[$i]\n";
               }
               $alignseqR2=$alignseqR2."$readseqR2[$i]";
            }
         }
         $alignseqR2=$padN.$alignseqR2;

         print STDERR "\t$countread" if ($countread%1000==0);
         print LOGFILE "\t$countread" if ($countread%1000==0);
         my $tostoreR1="$readidR1\tR1\t$nummatchR1\t$nummismatchR1\t$numctR1\t$numgaR1\t".($numothersR1+$numlqR1)."\t$readseqR1\t$atpositionsR1\t$agpositionsR1\t$acpositionsR1\t$tapositionsR1\t$tgpositionsR1\t$tcpositionsR1\t$gapositionsR1\t$gtpositionsR1\t$gcpositionsR1\t$capositionsR1\t$ctpositionsR1\t$cgpositionsR1\t$lqpositionsR1";
         my $tostoreR2="$readidR2\tR2\t$nummatchR2\t$nummismatchR2\t$numctR2\t$numgaR2\t".($numothersR2+$numlqR2)."\t$readseqR2\t$atpositionsR2\t$agpositionsR2\t$acpositionsR2\t$tapositionsR2\t$tgpositionsR2\t$tcpositionsR2\t$gapositionsR2\t$gtpositionsR2\t$gcpositionsR2\t$capositionsR2\t$ctpositionsR2\t$cgpositionsR2\t$lqpositionsR2";

	 # Condition cutoff rules below
         my $cond1pass="Cond1Fail";
         ($cond1pass="Cond1Pass") if (($nummatchR1>=$R1matchcut)&&($nummatchR2>=$R2matchcut));
         my $cond2pass="Cond2Fail";
         ($cond2pass="Cond2Pass") if ((($numctR1>$R1ctcutA)&&($numctR2>$R2ctcutA))&&(($numgaR1<=$R1gacutA)&&($numgaR2<=$R2gacutA))&&(($numothersR1+$numlqR1<=$R1otherscut)&&($numothersR2+$numlqR2<=$R2otherscut)));
         ($cond2pass="Cond2Pass") if ((($numgaR1>$R1gacutB)&&($numgaR2>$R2gacutB))&&(($numctR1<=$R1ctcutB)&&($numctR2<=$R2ctcutB))&&(($numothersR1+$numlqR1<=$R1otherscut)&&($numothersR2+$numlqR2<=$R2otherscut)));

         my $cat="$cond1pass\_$cond2pass";

	 print OUTFILE "$countpairs\t$tostoreR1\t$cat\n$countpairs\t$tostoreR2\t$cat\n";

         my $tostore3=">>>$countread\t>>>$readidR2\t>>>$cat\n";
	 $tostore3=$tostore3."R1: $atpositionsR1 $agpositionsR1 $acpositionsR1\n";
	 $tostore3=$tostore3."R1: $tapositionsR1 $tgpositionsR1 $tcpositionsR1\n";
	 $tostore3=$tostore3."R1: $gapositionsR1 $gtpositionsR1 $gcpositionsR1\n";
	 $tostore3=$tostore3."R1: $capositionsR1 $ctpositionsR1 $cgpositionsR1\nR1lq: $lqpositionsR1\n";
	 $tostore3=$tostore3."R2: $atpositionsR2 $agpositionsR2 $acpositionsR2\n";
	 $tostore3=$tostore3."R2: $tapositionsR2 $tgpositionsR2 $tcpositionsR2\n";
	 $tostore3=$tostore3."R2: $gapositionsR2 $gtpositionsR2 $gcpositionsR2\n";
	 $tostore3=$tostore3."R2: $capositionsR2 $ctpositionsR2 $cgpositionsR2\nR2lq: $lqpositionsR2\n\n";
         
	 for (my $j=0; $j<$reflen; ) {
            my $tempstart=$j+1;
            my $tempend=$j+length(substr($refseq,$j,50));
            $tostore3=$tostore3."Ref:      $tempstart\t".substr($refseq,$j,50)."\t$tempend\n";
            $tostore3=$tostore3."R1seq:    $tempstart\t".substr($alignseqR1,$j,50)."\t$tempend\n";
            $tostore3=$tostore3."R1qual:   $tempstart\t".substr($qualR1,$j,50)."\t$tempend\n";
            $tostore3=$tostore3."R2seq:    $tempstart\t".substr($alignseqR2,$j,50)."\t$tempend\n";
            $tostore3=$tostore3."R2qual:   $tempstart\t".substr($qualR2,$j,50)."\t$tempend\n\n";
            $j=$j+50
         }

         if ($cat eq "Cond1Pass_Cond2Pass") {
            $countPP++;
            print PP "$tostore3";
         }elsif ($cat eq "Cond1Pass_Cond2Fail") {
            $countPF++;
            print PF "$tostore3";
         }elsif ($cat eq "Cond1Fail_Cond2Pass") {
            $countFP++;
            print FP "$tostore3";
         }elsif ($cat eq "Cond1Fail_Cond2Fail") {
            $countFF++;
            print FF "$tostore3";
         }else {
            die "unknown category $cat...\n";
         }

      #re-initialize at end of loop
      $readseqR1=""; $readseqR2="";
      $nummatchR1=0; $nummatchR2=0;
      $nummismatchR1=0; $nummismatchR2=0;
      $numatR1=0; $numatR2=0;
      $numagR1=0; $numagR2=0;
      $numacR1=0; $numacR2=0;
      $numtaR1=0; $numtaR2=0;
      $numtgR1=0; $numtgR2=0;
      $numtcR1=0; $numtcR2=0;
      $numgaR1=0; $numgaR2=0;
      $numgtR1=0; $numgtR2=0;
      $numgcR1=0; $numgcR2=0;
      $numcaR1=0; $numcaR2=0;
      $numctR1=0; $numctR2=0;
      $numcgR1=0; $numcgR2=0;
      $numothersR1=0; $numothersR2=0;
      $numlqR1=0; $numlqR2=0;
      $atpositionsR1="ATpos"; $atpositionsR2="ATpos";
      $agpositionsR1="AGpos"; $agpositionsR2="AGpos";
      $acpositionsR1="ACpos"; $acpositionsR2="ACpos";
      $tapositionsR1="TApos"; $tapositionsR2="TApos";
      $tgpositionsR1="TGpos"; $tgpositionsR2="TGpos";
      $tcpositionsR1="TCpos"; $tcpositionsR2="TCpos";
      $gapositionsR1="GApos"; $gapositionsR2="GApos";
      $gtpositionsR1="GTpos"; $gtpositionsR2="GTpos";
      $gcpositionsR1="GCpos"; $gcpositionsR2="GCpos";
      $capositionsR1="CApos"; $capositionsR2="CApos";
      $ctpositionsR1="CTpos"; $ctpositionsR2="CTpos";
      $cgpositionsR1="CGpos"; $cgpositionsR2="CGpos";
      $otpositionsR1="OTpos"; $otpositionsR2="OTpos";
      $lqpositionsR1="LOWQUALpos"; $lqpositionsR2="LOWQUALpos";
      $alignseqR1=""; $alignseqR2="";
      }
   }
   print STDERR "\n$countread readpairs parsed...\n";
   print STDERR "$countpairs readpairs processed...\n";
   print LOGFILE "\n$countread readpairs parsed...\n";
   print LOGFILE "$countpairs readpairs processed...\n";
   close (OUTFILE);
   close (PP);
   close (PF);
   close (FP);
   close (FF);

   my $catsummaryout=$outputfile2;
   $catsummaryout=~s/\.txt$/\_Summary\.txt/gm;
   die "Cannot open $catsummaryout" unless (open(CATSUM, ">$catsummaryout"));
   print CATSUM "Category\tNumber\n";
   print CATSUM "Cond1Pass_Cond2Pass\t$countPP\n";
   print CATSUM "Cond1Pass_Cond2Fail\t$countPF\n";
   print CATSUM "Cond1Fail_Cond2Pass\t$countFP\n";
   print CATSUM "Cond1Fail_Cond2Fail\t$countFF\n";
   close(CATSUM);

   ### Tabulate Positional Results ###
   print STDERR "Part3 - Tabulating results\n";
   print LOGFILE "Part3 - Tabulating results\n";

   my $countline=0;
   my $countPP=0;
   my %combseqhash=();
   my %mismatchhash=();
   open (INPUT,"$outputfile1") or die "Cannot open xxx.01_allreads_mismatch_stats.txt file!!! \n";
   while (my $line=<INPUT>) {
      chomp $line;
      my $line2=$line;
      my $line3=$line;

      my @array=split(/\s+/,$line2);
      my $arraylen=scalar(@array);
      die "array len error...$arraylen...\n" if ($arraylen != 23);
      next if ($array[0] eq "number");
      $countline++;

      my $readnum=$array[0];
      my $readid=$array[1];
      my $readorientation=$array[2];
      my $readnummatch=$array[3];
      my $readnummismatch=$array[4];
      my $readnumct=$array[5];
      my $readnumga=$array[6];
      my $readnumothers=$array[7];
      my $readseq=$array[8];
      my $readatpositions=$array[9];
      my $readagpositions=$array[10];
      my $readacpositions=$array[11];
      my $readtapositions=$array[12];
      my $readtgpositions=$array[13];
      my $readtcpositions=$array[14];
      my $readgapositions=$array[15];
      my $readgtpositions=$array[16];
      my $readgcpositions=$array[17];
      my $readcapositions=$array[18];
      my $readctpositions=$array[19];
      my $readcgpositions=$array[20];
      my $readlqpositions=$array[21];
      my $readcat=$array[22];

      if ($readcat eq "Cond1Pass_Cond2Pass") {
         $countPP++;
         if ($readorientation eq "R1") {
            die "ERROR1!\n" if (exists $combseqhash{$readid});
            $combseqhash{$readid}=$readseq;
            $mismatchhash{"$readid\_R1"}="$readatpositions\:\:$readagpositions\:\:$readacpositions\:\:$readtapositions\:\:$readtgpositions\:\:$readtcpositions\:\:$readgapositions\:\:$readgtpositions\:\:$readgcpositions\:\:$readcapositions\:\:$readctpositions\:\:$readcgpositions\:\:OTpos\:\:$readlqpositions\:\:";
         }elsif ($readorientation eq "R2") {
            die "ERROR2!\n" if (!exists $combseqhash{$readid});
            $combseqhash{$readid}="$combseqhash{$readid}\/\/$readseq";
            $mismatchhash{"$readid\_R2"}="$readatpositions\:\:$readagpositions\:\:$readacpositions\:\:$readtapositions\:\:$readtgpositions\:\:$readtcpositions\:\:$readgapositions\:\:$readgtpositions\:\:$readgcpositions\:\:$readcapositions\:\:$readctpositions\:\:$readcgpositions\:\:OTpos\:\:$readlqpositions\:\:";
         }else {
            die "unknown orientation...$readorientation...\n";
         }
      }
      #last if ($countPP==400);
   }
   print STDERR "$countline lines parsed... $countPP PP...\n";
   print LOGFILE "$countline lines parsed... $countPP PP...\n";
   #sleep(3);

   my $countread=0;
   my %deduphash=();
   foreach my $readid (sort {$a cmp $b} (keys %combseqhash)) {
      $countread++;
      my $combseq=$combseqhash{$readid};
      if (!exists $deduphash{$combseq}) {
         $deduphash{$combseq}=$readid;
      }else {
         $deduphash{$combseq}="$deduphash{$combseq}\/\/$readid";
      }
   }
   my $numuniqseqs=scalar(keys(%deduphash));
   print STDERR "$countread PP readseqs parsed, $numuniqseqs are unique (if deduped)\n";
   print LOGFILE "$countread PP readseqs parsed, $numuniqseqs are unique (if deduped)\n";



   #initialize arrays for output
   my @Position=("Position");
   my @Nucleotide=("Nucleotide");
   my @CT_NumSeqs=("CTgroup_Num_Seqs");
   my @CT_NumA=("CTgroup_Num_A"); my @CT_NumT=("CTgroup_Num_T");
   my @CT_NumG=("CTgroup_Num_G"); my @CT_NumC=("CTgroup_Num_C");
   my @CT_NumOT=("CTgroup_Num_OthersNLowqual");
   my @CT_NumCT=("CTgroup_Num_CT"); my @CT_NumGA=("CTgroup_Num_GA");
   my @GA_NumSeqs=("GAgroup_Num_Seqs");
   my @GA_NumA=("gagroup_Num_A"); my @GA_NumT=("GAgroup_Num_T");
   my @GA_NumG=("GAgroup_Num_G"); my @GA_NumC=("GAgroup_Num_C");
   my @GA_NumOT=("GAgroup_Num_OthersNLowqual");
   my @GA_NumCT=("GAgroup_Num_CT"); my @GA_NumGA=("GAgroup_Num_GA");
   my @CTGA_NumSeqs=("CTGAgroup_Num_Seqs");
   my @CTGA_NumA=("CTGAgroup_Num_A"); my @CTGA_NumT=("CTGAgroup_Num_T");
   my @CTGA_NumG=("CTGAgroup_Num_G"); my @CTGA_NumC=("CTGAgroup_Num_C");
   my @CTGA_NumOT=("CTGAgroup_Num_OthersNLowqual");
   my @CTGA_NumCT=("CTGAgroup_Num_CT"); my @CTGA_NumGA=("CTGAgroup_Num_GA");
   my @ALL_NumSeqs=("ALL_Num_Seqs");
   my @ALL_NumA=("ALL_Num_A"); my @ALL_NumT=("ALL_Num_T");
   my @ALL_NumG=("ALL_Num_G"); my @ALL_NumC=("ALL_Num_C");
   my @ALL_NumOT=("ALL_Num_OthersNLowqual");
   my @ALL_NumCT=("ALL_Num_CT"); my @ALL_NumGA=("ALL_Num_GA");
   my @REF_NumA=("REF_Num_A"); my @REF_NumT=("REF_Num_T");
   my @REF_NumG=("REF_Num_G"); my @REF_NumC=("REF_Num_C");
   my @REF_NumOT=("REF_Num_OthersNLowqual");
   my @REF_NumCT=("REF_Num_CT"); my @REF_NumGA=("REF_Num_GA");
   foreach (my $i=1;$i<=$reflen;$i++) {
      push(@Position,$i);
      push(@Nucleotide,$refseq[$i-1]);
      push(@CT_NumA,0); push(@CT_NumT,0);
      push(@CT_NumG,0); push(@CT_NumC,0);
      push(@CT_NumOT,0);
      push(@CT_NumCT,0); push(@CT_NumGA,0);
      push(@GA_NumA,0); push(@GA_NumT,0);
      push(@GA_NumG,0); push(@GA_NumC,0);
      push(@GA_NumOT,0);
      push(@GA_NumCT,0); push(@GA_NumGA,0);
      push(@CTGA_NumA,0); push(@CTGA_NumT,0);
      push(@CTGA_NumG,0); push(@CTGA_NumC,0);
      push(@CTGA_NumOT,0);
      push(@CTGA_NumCT,0); push(@CTGA_NumGA,0);
      push(@ALL_NumA,0); push(@ALL_NumT,0);
      push(@ALL_NumG,0); push(@ALL_NumC,0);
      push(@ALL_NumOT,0);
      push(@ALL_NumCT,0); push(@ALL_NumGA,0);
      push(@REF_NumOT,0);
      push(@REF_NumCT,0); push(@REF_NumGA,0);
      if ($refseq[$i-1] eq "A") {
         push(@REF_NumA,1); push(@REF_NumT,0);
         push(@REF_NumG,0); push(@REF_NumC,0);
      }elsif ($refseq[$i-1] eq "T") {
         push(@REF_NumA,0); push(@REF_NumT,1);
         push(@REF_NumG,0); push(@REF_NumC,0);
      }elsif ($refseq[$i-1] eq "G") {
         push(@REF_NumA,0); push(@REF_NumT,0);
         push(@REF_NumG,1); push(@REF_NumC,0);
      }elsif ($refseq[$i-1] eq "C") {
         push(@REF_NumA,0); push(@REF_NumT,0);
         push(@REF_NumG,0); push(@REF_NumC,1);
      }else {
         die "ERROR5!\n";
      }
   }

   my $outputfile3="$sample\.03_Cond1Pass_Cond2Pass\.DuplicateCount\.txt";
   unless (open(OUTFILE1, ">$outputfile3")){
      print STDERR "Cannot open file \"$outputfile3\" to write to!!\n\n";
      exit;
   }
   print OUTFILE1 "R1_R2_Sequence\tNumCopies\tGroup\n";
   my $countdedup=0;
   my $countseq=0;
   my $countCTgroup=0;
   my $countGAgroup=0;
   my $countCTGAgroup=0;
   foreach my $key (sort {$a cmp $b} (keys %deduphash)) {
      $countdedup++;
      my $readid=$deduphash{$key};
      my $numdup=scalar(split(/\/\//,$readid));
      $readid=~s/\/\/.+//gm;

      #print STDERR "\n>>>$countdedup\t$readid...\tnumdups=$numdup\n";
      die "$readid\_R1 not found!!!\n" if (!exists $mismatchhash{"$readid\_R1"});
      die "$readid\_R2 not found!!!\n" if (!exists $mismatchhash{"$readid\_R2"});
      #print STDERR $mismatchhash{"$readid\_R1"}."\n";
      #print STDERR $mismatchhash{"$readid\_R2"}."\n";
      #print STDERR "$CT_R1  $GA_R1  $OT_R1  $LQ_R1\n";
      #print STDERR "$CT_R2  $GA_R2  $OT_R2  $LQ_R2\n";

      my %temphash=();
      #print STDERR "$readid\n";
      %temphash=fillhash(\%temphash,$mismatchhash{"$readid\_R1"},$mismatchhash{"$readid\_R2"});
      delete $temphash{''};
      #print STDERR "\n";

      my @Temp_NumA=@REF_NumA; my @Temp_NumT=@REF_NumT;
      my @Temp_NumG=@REF_NumG; my @Temp_NumC=@REF_NumC;
      my @Temp_NumOT=@REF_NumOT;
      my @Temp_NumCT=@REF_NumCT; my @Temp_NumGA=@REF_NumGA;
      #print STDERR "@Temp_NumA\n@Temp_NumT\n@Temp_NumG\n@Temp_NumC\n@Temp_NumOT\n\n";

      my $tempnumCT=0;
      my $tempnumGA=0;
      foreach my $pos (sort {$a <=> $b} keys(%temphash)) {
         my $mut=$temphash{$pos};
	 my ($mut1,$mut2)=split(//,$mut);

         if ($mut1 eq "A") {
            $Temp_NumA[$pos]=$Temp_NumA[$pos]-1;
         }elsif ($mut1 eq "T") {
            $Temp_NumT[$pos]=$Temp_NumT[$pos]-1;
         }elsif ($mut1 eq "G") {
            $Temp_NumG[$pos]=$Temp_NumG[$pos]-1;
         }elsif ($mut1 eq "C") {
            $Temp_NumC[$pos]=$Temp_NumC[$pos]-1;
         }elsif (($mut1 eq "L")||($mut1 eq "M")) {
            $Temp_NumA[$pos]=0 if ($Temp_NumA[$pos]==1);
            $Temp_NumT[$pos]=0 if ($Temp_NumT[$pos]==1);
            $Temp_NumG[$pos]=0 if ($Temp_NumG[$pos]==1);
            $Temp_NumC[$pos]=0 if ($Temp_NumC[$pos]==1);
	 }else {
            die "ERROR6! $mut1 $mut2\n";
	 }
         if ($mut2 eq "A") {
            $Temp_NumA[$pos]=$Temp_NumA[$pos]+1;
         }elsif ($mut2 eq "T") {
            $Temp_NumT[$pos]=$Temp_NumT[$pos]+1;
         }elsif ($mut2 eq "G") {
            $Temp_NumG[$pos]=$Temp_NumG[$pos]+1;
         }elsif ($mut2 eq "C") {
            $Temp_NumC[$pos]=$Temp_NumC[$pos]+1;
         }elsif (($mut2 eq "Q")||($mut2 eq "M")) {
            $Temp_NumOT[$pos]=$Temp_NumOT[$pos]+1;
	 }else {
            die "ERROR7! $mut1 $mut2\n";
	 }

	 if ($mut eq "CT") {
	    $tempnumCT++;
	    $Temp_NumCT[$pos]=$Temp_NumCT[$pos]+1;
         }elsif ($mut eq "GA") {
	    $tempnumGA++;
	    $Temp_NumGA[$pos]=$Temp_NumGA[$pos]+1;
         }
	 #print STDERR "   $pos $mut $mut1 $mut2 $tempnumCT $tempnumGA\n";
      }
      #print STDERR "@Temp_NumA\n@Temp_NumT\n@Temp_NumG\n@Temp_NumC\n@Temp_NumOT\n\n";

      my $increment=1;
      $increment=$numdup if ($removedup eq "No");

      $countseq=$countseq+$increment;
      foreach (my $i=1;$i<=$reflen;$i++) {
         $ALL_NumA["$i"]=$ALL_NumA["$i"]+($Temp_NumA["$i"]*$increment);
         $ALL_NumT["$i"]=$ALL_NumT["$i"]+($Temp_NumT["$i"]*$increment);
         $ALL_NumG["$i"]=$ALL_NumG["$i"]+($Temp_NumG["$i"]*$increment);
         $ALL_NumC["$i"]=$ALL_NumC["$i"]+($Temp_NumC["$i"]*$increment);
         $ALL_NumOT["$i"]=$ALL_NumOT["$i"]+($Temp_NumOT["$i"]*$increment);
         $ALL_NumCT["$i"]=$ALL_NumCT["$i"]+($Temp_NumCT["$i"]*$increment);
         $ALL_NumGA["$i"]=$ALL_NumGA["$i"]+($Temp_NumGA["$i"]*$increment);
      }

      if ($tempnumCT > $tempnumGA) {
         # we define this as the CT subgroup seqs
         $countCTgroup=$countCTgroup+$increment;
         foreach (my $i=1;$i<=$reflen;$i++) {
            $CT_NumA["$i"]=$CT_NumA["$i"]+($Temp_NumA["$i"]*$increment);
            $CT_NumT["$i"]=$CT_NumT["$i"]+($Temp_NumT["$i"]*$increment);
            $CT_NumG["$i"]=$CT_NumG["$i"]+($Temp_NumG["$i"]*$increment);
            $CT_NumC["$i"]=$CT_NumC["$i"]+($Temp_NumC["$i"]*$increment);
            $CT_NumOT["$i"]=$CT_NumOT["$i"]+($Temp_NumOT["$i"]*$increment);
            $CT_NumCT["$i"]=$CT_NumCT["$i"]+($Temp_NumCT["$i"]*$increment);
            $CT_NumGA["$i"]=$CT_NumGA["$i"]+($Temp_NumGA["$i"]*$increment);
         }
         print OUTFILE1 "$key\t$numdup\tCTgroup\n";
      }elsif ($tempnumCT < $tempnumGA) {
         # we define this as the GA subgroup seqs
         $countGAgroup=$countGAgroup+$increment;
         foreach (my $i=1;$i<=$reflen;$i++) {
            $GA_NumA["$i"]=$GA_NumA["$i"]+($Temp_NumA["$i"]*$increment);
            $GA_NumT["$i"]=$GA_NumT["$i"]+($Temp_NumT["$i"]*$increment);
            $GA_NumG["$i"]=$GA_NumG["$i"]+($Temp_NumG["$i"]*$increment);
            $GA_NumC["$i"]=$GA_NumC["$i"]+($Temp_NumC["$i"]*$increment);
            $GA_NumOT["$i"]=$GA_NumOT["$i"]+($Temp_NumOT["$i"]*$increment);
            $GA_NumCT["$i"]=$GA_NumCT["$i"]+($Temp_NumCT["$i"]*$increment);
            $GA_NumGA["$i"]=$GA_NumGA["$i"]+($Temp_NumGA["$i"]*$increment);
         }
         print OUTFILE1 "$key\t$numdup\tGAgroup\n";
      }elsif ($tempnumCT == $tempnumGA) {
         # we define this as the CTGAequal subgroup seqs
         $countCTGAgroup=$countCTGAgroup+$increment;
         foreach (my $i=1;$i<=$reflen;$i++) {
            $CTGA_NumA["$i"]=$CTGA_NumA["$i"]+($Temp_NumA["$i"]*$increment);
            $CTGA_NumT["$i"]=$CTGA_NumT["$i"]+($Temp_NumT["$i"]*$increment);
            $CTGA_NumG["$i"]=$CTGA_NumG["$i"]+($Temp_NumG["$i"]*$increment);
            $CTGA_NumC["$i"]=$CTGA_NumC["$i"]+($Temp_NumC["$i"]*$increment);
            $CTGA_NumOT["$i"]=$CTGA_NumOT["$i"]+($Temp_NumOT["$i"]*$increment);
            $CTGA_NumCT["$i"]=$CTGA_NumCT["$i"]+($Temp_NumCT["$i"]*$increment);
            $CTGA_NumGA["$i"]=$CTGA_NumGA["$i"]+($Temp_NumGA["$i"]*$increment);
         }
         print OUTFILE1 "$key\t$numdup\tCTGAgroup\n";
      }else {
	 die "ERROR8!!!\n";
      }
   }
   foreach (my $i=1;$i<=$reflen;$i++) {
      push(@CT_NumSeqs,$countCTgroup);
      push(@GA_NumSeqs,$countGAgroup);
      push(@CTGA_NumSeqs,$countCTGAgroup);
      push(@ALL_NumSeqs,($countCTgroup+$countGAgroup+$countCTGAgroup));
   }
   print STDERR "$countCTgroup $countGAgroup $countCTGAgroup\n";

   my $outputfile4="$sample\.03_Cond1Pass_Cond2Pass\.PositionMutationCount\.txt";
   unless (open(OUTFILE2, ">$outputfile4")){
      print STDERR "Cannot open file \"$outputfile4\" to write to!!\n\n";
      exit;
   }
   
   print OUTFILE2 "@Position\n@Nucleotide\n\n";
   print OUTFILE2 "@CT_NumSeqs\n@CT_NumA\n@CT_NumT\n@CT_NumG\n@CT_NumC\n@CT_NumCT\n@CT_NumGA\n@CT_NumOT\n\n";
   print OUTFILE2 "@GA_NumSeqs\n@GA_NumA\n@GA_NumT\n@GA_NumG\n@GA_NumC\n@GA_NumCT\n@GA_NumGA\n@GA_NumOT\n\n";
   print OUTFILE2 "@CTGA_NumSeqs\n@CTGA_NumA\n@CTGA_NumT\n@CTGA_NumG\n@CTGA_NumC\n@CTGA_NumCT\n@CTGA_NumGA\n@CTGA_NumOT\n\n";
   print OUTFILE2 "@ALL_NumSeqs\n@ALL_NumA\n@ALL_NumT\n@ALL_NumG\n@ALL_NumC\n@ALL_NumCT\n@ALL_NumGA\n@ALL_NumOT\n\n";

   print STDERR "Positional mutation counts tabulated and printed...\n";
   print LOGFILE "Positional mutation counts tabulated and printed...\n";

   close(OUTFILE1);
   close(OUTFILE2);
   close (LOGFILE);
}
print STDERR ">>> Sequence Analysis completed...\n";

#TODO - Plotting
#   system("echo \"Part3 Plotting\" >> $logfile");
#   my $sample=$infile;
#   $sample=~s/^\S+\///gm;
#   $sample=~s/\.trimmed.fastq//gm;
#   my $cmd3="R < $script04 --no-save $outputdir $sample $hotstring $plotxmin $plotxmax $plotymin $plotymax >> $logfile 2>&1";
#   system("echo $cmd3 >> $logfile");
#   system($cmd3)==0 or die "\nFAILED TO EXECUTE: $cmd3\n";





