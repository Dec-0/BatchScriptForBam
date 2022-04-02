# Package Name
package Parameter::BinList;

# Exported name
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(BinListGet BinSearch DayStamp TimeStamp DayTimeStamp ScriptBegin TimeString Sec2String StorageAddComma RandString IfDirExist IfFileExist FileLock FileUnLock);

use FindBin qw($Bin);
use File::Basename;

sub BinListGet
{
	my $BinList = $Bin . "/.BinList.xls";
	return $BinList if(-s $BinList);
	
	$BinList = "/biobiggen/data/headQuarter/user/xiezhangdong/Scripts/.BinList.xls";
	die "[ Error ] File not exist ($BinList).\n" unless(-s $BinList);
	
	return $BinList;
}

sub BinSearch
{
	my ($BinName,$BinList,$JumpConfirmFlag) = @_;
	
	die "[ Error ] File not exist ($BinList).\n" unless(-s $BinList);
	my $Return = `grep '^$BinName'\$'\\t' $BinList`;
	chomp $Return;
	
	die "Name ($BinName) not unique.\n($Return)\n" if($Return =~ /\n/);
	#print "[ Info ] $Return\n";
	my @Cols = split /\t/, $Return;
	# 替换成Bin目录
	if($Cols[1] =~ s/^BINDIR//)
	{
		$Cols[1] = $Bin . $Cols[1];
	}
	# 替换成Scripts目录
	elsif($Cols[1] =~ s/^DEFAULTDIR//)
	{
		my @Items = split /\//, $Bin;
		my $LastId = -1;
		for($i = $#Items;$i >= 0;$i --)
		{
			if($Items[$i] eq "Scripts")
			{
				$LastId = $i - 1;
				last;
			}
			elsif($Items[$i] eq "xiezhangdong")
			{
				$LastId = $i;
				last;
			}
		}
		
		if($LastId >= 0)
		{
			my $RDir = join("/",@Items[0 .. $LastId]);
			$Cols[1] = $RDir . $Cols[1];
		}
	}
	die "[ Error ] $Cols[1] for $BinName not exist!\n" unless($Cols[1] && (-e $Cols[1] || $JumpConfirmFlag ));
	
	return $Cols[1];
}

sub DayStamp
{
	my @temp_time = localtime();
	my $localtime_year = $temp_time[5] + 1900;
	my $localtime_month = $temp_time[4] + 1;
	my $DayStamp = $localtime_year . "/" . $localtime_month . "/" . $temp_time[3];
	
	return $DayStamp;
}

sub TimeStamp
{
	my @temp_time = localtime();
	my $TimeStamp = $temp_time[2] . ":" . $temp_time[1] . ":" . $temp_time[0];
	
	return $TimeStamp;
}

sub DayTimeStamp
{
	my $DayStamp = &DayStamp();
	my $TimeStamp = &TimeStamp();
	my $Stamp = $DayStamp . " " . $TimeStamp;
	
	return $Stamp;
}

sub ScriptBegin
{
	my ($MuteFlag,$Name) = @_;
	
	$BeginTime = time;
	if($BeginTime && !$MuteFlag)
	{
		my $Stamp = &DayTimeStamp();
		print "[ $Stamp ] This script begins.\n" unless($Name);
		print "[ $Stamp ] This script ($Name) begins.\n" if($Name);
	}
	elsif(!$BeginTime)
	{
		die "[ Error ] Fail to get system time with 'time'.\n";
	}
	
	return $BeginTime;
}

sub TimeString
{
	my ($TimeCurrent,$TimeBegin) = @_;
	my $TimeString = "";
	
	my $Sec = $TimeCurrent - $TimeBegin;
	
	my $Day = int($Sec / 86400);
	$TimeString .= $Day . "d" if($Day);
	
	my $Tail = $Sec % 86400;
	my $Hour = int($Tail / 3600);
	$TimeString .= $Hour . "h" if($Hour);
	
	$Tail = $Tail % 3600;
	my $Min = int($Tail / 60);
	$TimeString .= $Min . "min" if($Min);
	
	$Sec = $Tail % 60;
	$TimeString .= $Sec . "s" if($Sec);
	
	$TimeString = "0s" unless($TimeString);
	
	return $TimeString;
}

sub Sec2String
{
	my $Sec = $_[0];
	my $TimeString = "";
	
	my $Hour = int($Sec / 3600);
	$Hour = "0" . $Hour unless(length($Hour) >= 2);
	
	$Sec = $Sec % 3600;
	my $Min = int($Sec / 60);
	$Min = "0" . $Min unless(length($Min) >= 2);
	
	$Sec = $Sec % 60;
	$Sec = "0" . $Sec unless(length($Sec) >= 2);
	
	$TimeString = join(":",$Hour,$Min,$Sec);
	
	return $TimeString;
}

sub StorageAddComma
{
	my $Mem = $_[0];
	
	my $RoundNum = 0;
	my @Items = split //, $Mem;
	for(my $i = $#Items;$i > 0;$i --)
	{
		$RoundNum ++;
		if($RoundNum == 3)
		{
			$Items[$i] = "," . $Items[$i];
			$RoundNum = 0;
		}
	}
	my $ConvertedMem = join("",@Items);
	
	return $ConvertedMem;
}

sub RandString
{
	my $Len = $_[0];
	
	die "[ Error ] The length of rand string is not pure number ($Len).\n" if($Len =~ /\D/);
	my @Items = ();
	for my $i (0 .. 9)
	{
		push @Items, $i;
	}
	for my $i (65 .. 90)
	{
		push @Items, chr($i);
	}
	for my $i (97 .. 122)
	{
		push @Items, chr($i);
	}
	my $Full = @Items;
	my $String = "";
	for my $i (1 .. $Len)
	{
		$String .= $Items[int(rand($Full))];
	}
	
	return $String;
}

sub IfDirExist
{
	my @Dir = @_;
	
	for my $i (0 .. $#Dir)
	{
		$Dir[$i] =~ s/\/$//;
		die "[ Error ] Directory ($Dir[$i]) not exist.\n" unless($Dir[$i] && -d $Dir[$i]);
	}
	
	if($#Dir > 0)
	{
		return @Dir;
	}
	elsif($#Dir == 0)
	{
		return $Dir[0];
	}
	else
	{
		return "";
	}
}

sub IfFileExist
{
	my @File = @_;
	
	for my $i (0 .. $#File)
	{
		die "[ Error ] File ($File[$i]) not exist.\n" unless($File[$i] && -e $File[$i]);
	}
	
	return 1;
}

sub FileLock
{
	my $File = $_[0];
	my $SleepTime = $_[1];
	$SleepTime = 60 unless($SleepTime);
	
	my $Dir = dirname $File;
	die "[ Error ] Directory($Dir) not exist.\n" unless(-d $Dir);
	my $BaseName = basename $File;
	my $LockFile = $File . ".lock";
	
	my $RandString = "";
	while(1)
	{
		# 需要其它程序删除这个lock文件;
		unless(-e $LockFile)
		{
			# 新建lock文件;
			$RandString = &RandString(30);
			`echo '$RandString' > $LockFile`;
			
			# 确认该lock文件是我们新建的，而不是其它程序的;
			my $Return = `cat $LockFile`;
			chomp $Return;
			last if($Return eq $RandString);
		}
		
		#print "[ Info ] ( $BaseName ) File locked, waiting ...\n";
		sleep($SleepTime);
	}
	
	return $RandString;
}

sub FileUnLock
{
	my $File = $_[0];
	my $RandString = $_[1];
	
	my $LockFile = $File . ".lock";
	die "[ Error ] Lock file ($LockFile) not exist.\n" unless(-s $LockFile);
	my $Return = `cat $LockFile`;
	chomp $Return;
	if($Return eq $RandString)
	{
		`rm $LockFile`;
	}
	else
	{
		die "[ Error ] Content in Lock file not consistent with $RandString.\n";
	}
	
	return 1;
}

1;