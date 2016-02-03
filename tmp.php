<?php


$subdir = rand();
$cmd = 'mkdir output/' . $subdir;
shell_exec ( $cmd );
$refdir = 'ref/';

$command="./exomefilter -numRefIndv 2535 -genoFile $refdir/chris.reduced.genos -chromFile $refdir/chris.reduced.chrom -posFile $refdir/chris.reduced.pos -out " . $subdi


if( $_POST["Design"]=="unrelated"){$command = $command . " -unrelatedStudy " . $_POST["NumUnrelated"];}

//echo $_POST["pedFileData"];
if( $_POST["Design"]=="family") {
        //echo $_POST["pedFileData"];
        //$target_path = "uploads/";$target_path .= basename( $_FILES['pedfile']['name']);
        //if(move_uploaded_file($_FILES['pedfile']['tmp_name'], $target_path)) {
        //      echo "upload successful";
        //} else {
        //      echo "There was an error uploading the file, please try again.";
        //}
        $file = "uploads/" . $subdir.".ped";
        //$current = file_get_contents($file);
        $current = $_POST["pedFileData"];
     $e =explode("\n",$current);

        for ($i=0;$i<count($e);$i++){$data .= rtrim($e[$i],"\r") . "\n";}
        //$current = rtrim($current);
        file_put_contents($file, $data);
        //echo $file;
}

$command .= " -pedigree " . $file;

$command .= " -" . $_POST["Filter1"] . " " . $_POST["Text1"];
$command .= " -" . $_POST["Filter2"] . " " . $_POST["Text2"];
$command .= " -" . $_POST["Filter3"] . " " . $_POST["Text3"];
$command .= " -" . $_POST["Filter4"] . " " . $_POST["Text4"];
$command .= " -" . $_POST["Filter5"] . " " . $_POST["Text5"];
$command .= " -" . $_POST["Filter6"] . " " . $_POST["Text6"];
$command .= " -" . $_POST["Filter7"] . " " . $_POST["Text7"];
$command .= " -" . $_POST["Filter8"] . " " . $_POST["Text8"];
$command .= " -" . $_POST["Filter9"] . " " . $_POST["Text9"];
$command .= " -" . $_POST["Filter10"] . " " . $_POST["Text10"];
$command .= " -" . $_POST["Filter11"] . " " . $_POST["Text11"];
$command .= " -" . $_POST["Filter12"] . " " . $_POST["Text12"];


$command .= " -rep " . $_POST["reps"];

shell_exec ( $command );


$rcmd = "Rscript output/" . $subdir . "/plots.R";

shell_exec ( $rcmd );

//header('Content-Type: application/json');
//$img = "output/" . $subdir . $_POST["Filter1"] . ".png";
//echo json_encode( $img );

$dirname = "output/" . $subdir . "/";
$images = glob($dirname."*.png");
echo '<br><br><br><br>Results<br><br><br>';
//foreach($images as $image) { echo '<img src="'.$image.'" /><br><br><br><br />'; }
for ($i=1;$i<count($images)+1;$i++){ $a = $dirname . $_POST["Filter$i"] . ".png"; echo '<img src="'.$a.'" /><br><br><br><br />'; }

?>
