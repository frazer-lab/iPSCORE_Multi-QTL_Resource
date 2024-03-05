#!/usr/bin/perl

use strict;

my $prv_name = "";
while(<STDIN>) {
	s/\s+$//;
	my($chr, $start, $end, $name) = (split /\t/)[0..3];
	$name =~ s/\D*$//;
	$name eq $prv_name && next;
	print join("\t", $name, $chr, $start, $end, "."), "\n";
	$prv_name = $name;
}

exit;
