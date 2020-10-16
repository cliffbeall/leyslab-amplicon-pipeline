<?php
ini_set('memory_limit', '12G');
//read fasta file and populate array
//fastaArray[seqNum][0]=>sequenceID[1]=>sequence
function explodeFastaFile()
{
	$numArgs=func_num_args();
	$arg=func_get_args($numArgs);
	
	$fastaFile=$arg[0];

	$fastaInput=fopen($fastaFile,"r");
	$fastaArray=array();
	$i=-1;

	while(!feof($fastaInput))
	{
		$line=fgets($fastaInput,8192);
		if(strncmp($line,">",1)==0)
		{
			$i++;
			$fastaArray[$i]=array(0=>trim(substr($line,1)),1=>"");
		}
		else
		{
			$fastaArray[$i][1].=trim($line);
		}
	}
	
	fclose($fastaInput);
	return $fastaArray;
}

//recursive merge sort
function mergeSort($hits,$compare)
{
	if(count($hits)<2)
	{
		return $hits;
	}
	else
	{
		$mid=(int)(count($hits)/2);
		
		$slice1=array_slice($hits,0,$mid);
		$slice2=array_slice($hits,$mid);
		
		$slice1=mergeSort($slice1,$compare);
		$slice2=mergeSort($slice2,$compare);
		
		$merged=mergeMerge($slice1,$slice2,$compare);
		
		return $merged;
	}
}

//returns merged slices
//compare is name of the comparator f(x)
function mergeMerge($slice1,$slice2,$compare)
{
	$result=array();
	$totalNum=count($slice1)+count($slice2);
	$i=0;
	$j=0;
	while(count($result)<$totalNum)
	{

		if(@$slice1[$i]==NULL)
		{
			$result[count($result)]=$slice2[$j];
			$j++;
		}
		else if(@$slice2[$j]==NULL)
		{
			$result[count($result)]=$slice1[$i];
			$i++;
		}
		else if(call_user_func($compare,@$slice1[$i],@$slice2[$j]))
		{
			$result[count($result)]=$slice1[$i];
			$i++;
		}
		else
		{
			$result[count($result)]=$slice2[$j];
			$j++;
		}
	}
	return $result;
}

//returns true if hit1<=hit2
function compareHits($hit1,$hit2)
{
	if($hit1[1]>=$hit2[1])
	{
		return true;
	}
	else
	{

		return false;
	}
}

function baseMatch($base1,$base2)
{
	$base1=strtoupper($base1);
	$base2=strtoupper($base2);
	$i=0;
	
	while($i++<2)
	{
		if($base1=='-'||$base2=='-'){return false;}
		if($base1==$base2){return true;}
		if($base1=='N'){return true;}
		if($base1=='A'){if($base2=='W'||$base2=='H'||$base2=='D'||$base2=='M'||$base2=='R'||$base2=='V'||$base2=='N'){return true;}else{return false;}}
		if($base1=='T'){if($base2=='W'||$base2=='Y'||$base2=='H'||$base2=='K'||$base2=='D'||$base2=='B'||$base2=='N'){return true;}else{return false;}}
		if($base1=='G'){if($base2=='K'||$base2=='D'||$base2=='B'||$base2=='N'||$base2=='R'||$base2=='S'||$base2=='V'){return true;}else{return false;}}
		if($base1=='C'){if($base2=='Y'||$base2=='H'||$base2=='B'||$base2=='N'||$base2=='M'||$base2=='S'||$base2=='V'){return true;}else{return false;}}
		if($base1=='V')
		{
			if($base2=='T'){return false;}else{return true;}
		}
		if($base1=='W')
		{
			if($base2=='S'||$base2=='G'||$base2=='C'){return false;}else{return true;}
		}
		if($base1=='Y')
		{
			if($base2=='R'||$base2=='A'||$base2=='G'){return false;}else{return true;}
		}
		if($base1=='H')
		{
			if($base2=='G'||$base2=='-'){return false;}else{return true;}
		}
		if($base1=='K') 
		{
			if($base2=='M'||$base2=='A'||$base2=='C'){return false;}else{return true;}
		}
		if($base1=='D')
		{
			if($base2=='C'){return false;}else{return true;}
		}
		if($base1=='B')
		{
			if($base2=='A'){return false;}else{return true;}
		}
		if($base1=='M')
		{
			if($base2=='K'||$base2=='T'||$base2=='G'){return false;}else{return true;}
		}
		if($base1=='R')
		{
			if($base2=='Y'||$base2=='T'||$base2=='C'){return false;}else{return true;}
		}
		if($base1=='S')
		{
			if($base2=='W'||$base2=='A'||$base2=='T'){return false;}else{return true;}
		}		
		$tmp=$base2;
		$base2=$base1;
		$base1=$tmp;
	}
	
	return false;
}

//Get % identity
//counting ambigious bases as match if they can
//
function alignmentStats($alignQuerySeq,$alignHitSeq,&$hitIdentity,&$hitAlignmentLength)
{
	$numMatch=0;
	$hitAlignmentLength=strlen($alignQuerySeq);

	for($i=0;$i<$hitAlignmentLength;$i++)
	{
		if(baseMatch($alignQuerySeq[$i],$alignHitSeq[$i]))
		{
			$numMatch++;
		}	
	}

	$hitIdentity=$numMatch/$hitAlignmentLength;

}

//reads tab result file and populates
//an array with the data
//resultArray[queryId][0]=>hitAcc.[1]=>hit%Ident.[2]=>alignLength
function populateResultArray($tabResultInput,$minAlignmentLength)
{
	$resultArray=array();
	
	while(!feof($tabResultInput))
	{
		$line=fgets($tabResultInput,8192);
		if($line=="")
		{
			break;
		}
		$line=explode("\t",$line);
		
		$queryId=$line[0];
		$hitAccession=$line[1];
		$alignQuerySeq=$line[2];
		$alignHitSeq=$line[3];

		alignmentStats($alignQuerySeq,$alignHitSeq,$hitIdentity,$hitAlignmentLength);

		if($hitAlignmentLength<$minAlignmentLength)
		{
			$hitIdentity/=10;
		}
		
		
		
		if(array_key_exists($queryId,$resultArray))
		{
			$newHitIndex=count($resultArray[$queryId]);
			$resultArray[$queryId][$newHitIndex]=array(0=>$hitAccession,1=>$hitIdentity,2=>$hitAlignmentLength);
		}
		else
		{
			$resultArray[$queryId]=array();
			$resultArray[$queryId][0]=array(0=>$hitAccession,1=>$hitIdentity,2=>$hitAlignmentLength);
		}
	}
	return $resultArray;
}

function buildTextOut($resultArray,$fastaArray,$maxNumHitsToReturn,$padResults,$showLength)
{
	$resultText="";
	
	foreach($fastaArray as $query)
	{
		$queryName=$query[0];
		$querySequence=$query[1];
		
		if(array_key_exists($queryName,$resultArray))
		{
			$numHits=count($resultArray[$queryName]);
			$resultText.=$queryName;
			
			for($i=0;$i<$maxNumHitsToReturn;$i++)
			{
				if($i<$numHits)
				{
					$hitName=$resultArray[$queryName][$i][0];
					$hitPercId=$resultArray[$queryName][$i][1];
					$hitLength=$resultArray[$queryName][$i][2];
					$resultText.="\t$hitName\t$hitPercId";
					if($showLength)
					{
						$resultText.="\t$hitLength";
					}
				}
				elseif($padResults)
				{
					$resultText.="\t\t";
					if($showLength)
					{
						$resultText.="\t";
					}
				}
			}
			$resultText.="\t$querySequence\n";
		}
		else
		{
			$resultText.=$queryName."\t[NO_HIT]\t0";
			if($padResults)
			{
				$resultText.=str_repeat("\t\t",$maxHitsToReturn-1);	
			}
			$resultText.="\t$querySequence\n";
		}
	}
	return $resultText;
}
function summarize()
{
	$numArgs=func_num_args();
	$arg=func_get_args($numArgs);
	
	$fastaFile=$arg[0];
	$tabResultFile=$arg[1];
	$maxNumHitsToReturn=$arg[2];
	$minAlignmentLength=$arg[3];
	$padResults=($numArgs>4)?$arg[4]:false;
	$padResults=(@$padResult=='t'||@$padResult=='T')?true:false;
	$showLength=($numArgs>5)?$arg[5]:false;
	$showLength=($showLength=='t'||$showLength=='T')?true:false;
	
	$fastaArray=explodeFastaFile($fastaFile);
	$tabResultInput=fopen($tabResultFile,"r");

	$resultArray=populateResultArray($tabResultInput,$minAlignmentLength);

	foreach($resultArray as &$queryHits)
	{
		$queryHits=mergeSort($queryHits,'compareHits');
	}
	
	$resultText=buildTextOut($resultArray,$fastaArray,$maxNumHitsToReturn,$padResults,$showLength);
	
	fclose($tabResultInput);
	
	return $resultText;
	
}
	
//php ncbiSummarize.php QUERY.FA RESULTS.TAB MAXHITS MINALIGNLENGTH PADRESULTS SHOWLENGTH OUTFILE
file_put_contents($argv[7],summarize($argv[1],$argv[2],$argv[3],$argv[4],$argv[5],$argv[6]));

//$testArray=array(0=>array(0=>"A",1=>2),1=>array(0=>"B",1=>3),2=>array(0=>"C",1=>3),3=>array(0=>"D",1=>4),4=>array(0=>"E",1=>1));
//$testArray=mergeSort($testArray,'compareHits');
//foreach($testArray as $hit)
//{
//	echo $hit[0]."-".$hit[1]."\n";
//}

//var_dump($testArray);
//$a=array('A','T','G','C','W','Y','H','K','D','B','M','R','S','V','N','-');
//$b=array('A','T','G','C','W','Y','H','K','D','B','M','R','S','V','N','-');
	
//for($i=0;$i<count($a);$i++)
//{
//	for($j=0;$j<count($b);$j++)
//	{
//		echo $a[$i]."=".$b[$j]."?".(baseMatch($a[$i],$b[$j])?"true\n":"false\n");
//	}
//}
?>
