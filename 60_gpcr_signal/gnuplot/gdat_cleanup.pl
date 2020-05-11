#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';

while (<>) {
    chomp;
    my @fields = split;
    if ($fields[0] eq "#") {
        $fields[0]="%";
        print join(" ",@fields);
        print("\n");
    } else {
        my @fields_str = map {sprintf("%.3e",$_)} @fields;
        print join(" ",@fields_str);
        print("\n");

    }


}


1;