#!/usr/bin/perl
#
#              INGL?S/ENGLISH
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  http://www.gnu.org/copyleft/gpl.html
#
#
#             PORTUGU?S/PORTUGUESE
#  Este programa ? distribu?do na expectativa de ser ?til aos seus
#  usu?rios, por?m N?O TEM NENHUMA GARANTIA, EXPL?CITAS OU IMPL?CITAS,
#  COMERCIAIS OU DE ATENDIMENTO A UMA DETERMINADA FINALIDADE.  Consulte
#  a Licen?a P?blica Geral GNU para maiores detalhes.
#  http://www.gnu.org/copyleft/gpl.html
#
#  Copyright (C) 2007  Funda??o Hemocentro de Ribeir?o Preto
#
#  Laborat?rio de Bioinform?tica
#  BiT -  Bioinformatic Team
#  Funda?o Hemocentro de Ribeir?o Preto
#  Rua Tenente Cat?o Roxo, 2501
#  Ribeir?o Preto - S?o Paulo
#  Brasil
#  CEP 14051-140
#  Fone: 55 16 39639300 Ramal 9603
#
#  Matheus B?rger
#  matheus@lgmb.fmrp.usp.br
#  http://lgmb.fmrp.usp.br
#
# $Id$

=head1 NAME 

=head1 SYNOPSIS

=head1 ABSTRACT

=head1 DESCRIPTION

=head1 AUTHOR

Matheus B?rger E<lt>matheus@lgmb.fmrp.usp.brE<gt>

Copyright (c) 2007 Regional Blood Center of Ribeir?o Preto

=head1 LICENSE

GNU General Public License

http://www.gnu.org/copyleft/gpl.html


=cut


use strict;
use warnings;
use Getopt::Long;
use Digest::MD5;

#Usage("Too few arguments") if $#ARGV<0;
GetOptions(	"h|?|help"=> sub { Usage(); } ) or Usage();
	   
#my $html = "tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/lgg/cgcc/jhu-usc.edu/";
#
my $base = "https://tcga-data.nci.nih.gov/";
#
#my $batchName = "humanmethylation450";

my $html = "tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/lgg/cgcc/jhu-usc.edu/humanmethylation450/methylation/";

`wget --no-check-certificate $base$html -O /tmp/teste.html`;


open(my $AUX, "<", "/tmp/teste.html");


my %files;

while(my $line = <$AUX>){
    chomp($line);
    if($line =~ /<a href=\"(\S+.tar.gz\S*)\">/){
        my $filename = $1;
        if($filename =~ /(\S+)\.Level_3\.(\d+)\./){
            #print $filename,"\n";
            my $key = $1."/".$2;
            $files{$key} = {} unless exists $files{$key};
            if($filename =~ /\.md5/){
                $files{$key}->{'md5'} = $filename;
            }elsif($filename =~ /\.tar\.gz/){
                $files{$key}->{'tgz'} = $filename;
            }
        }
    }
}

foreach my $k (keys %files){
    #print $k, "\n";
    my $count = 0;
    do{
        print `wget -c --no-check-certificate $base/$html/$files{$k}->{"tgz"}`;
        print `wget -c --no-check-certificate $base/$html/$files{$k}->{"md5"}`;
        my $ctx = Digest::MD5->new;
        open(my $NULL, "<", $files{$k}->{"tgz"});
        $ctx->addfile($NULL);
        $files{$k}->{'my_md5'} = $ctx->hexdigest;
        open(my $MD5, '<', $files{$k}->{"md5"});
        my $aux = <$MD5>;
        chomp($aux);
        my @md5 = split(/\s+/,$aux);
        close($MD5); 
        $files{$k}->{'their_md5'} = $md5[0];
        #print join("\t", @md5), "\n";
        #print $files{$k}->{'my_md5'},"\n";
        $count++;
        #print ">",$files{$k}->{'their_md5'},"<\t>",$files{$k}->{'my_md5'},"<\n";
    }while($files{$k}->{'their_md5'} ne  $files{$k}->{'my_md5'} and $count <= 4);
    if($files{$k}->{'their_md5'} ne  $files{$k}->{'my_md5'}){
        $files{$k}->{'error'} = "YES";
    }else{
        $files{$k}->{'error'} = "NO";
    }
}

open(my $OUT, ">", "info.txt");
select($OUT);
foreach my $k (keys %files){
    print $k,"\t";
    foreach my $kk (keys %{$files{$k}}){
        print $kk,':',$files{$k}->{$kk},"\t";
    }
    print "\n";
}
close($OUT);
#my @files_s = sort @files;
#
##print join("\n", @files_s);
#
#my $d_file = pop(@files_s);
#
#print `wget --no-check-certificate $base/$html/$d_file.tar.gz`;




# Subroutines 

sub Usage {
	my($msg) = @_;
	print STDERR "\nERR: $msg\n\n" if $msg;
	print STDERR qq[$0  ].q[$Revision$].qq[\n];
	print STDERR<<EOU;
Matheus B?rger (matheus\@bit.fmrp.usp.br) 
(c)2007 Regional Blood Center of Ribeir?o Preto

Usage 

	$0	

Argument(s)

	-h	--help	Help

EOU
	exit(1);
}


