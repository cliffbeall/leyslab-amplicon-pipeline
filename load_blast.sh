# load_blast.sh {run} {blast_db}
# to be run in a folder containing:
#	1. output from blast.bsh in a folder named blast
#	2. an excel file with plate information in the format Run*.xlsx
#	3. a groups file with the format {run}.q28.groups	

if [ -z "$2" ]; then
	blast_db="data.core_vaginal"
	else
	blast_db=$2
	fi




#drop if exists
echo "CREATE DATABASE IF NOT EXISTS $1" | mysql -u root
echo "DROP TABLE IF EXISTS max_scores,groups,plate" | mysql -u root $1
#create tables from mysqldb bioschemas
echo "CREATE TABLE max_scores like bioschemas.max_scores; CREATE TABLE groups like bioschemas.groups; CREATE TABLE plate like bioschemas.plate" | mysql -u root $1



#link groups file
ln -s $1.q28.groups groups

#find Run excel file
plate_file=$(ls | grep Run | grep xlsx);


#format plate in R and get study from sample name
echo 'library(gdata);
	map <- read.xls("'$plate_file'",sheet=1);
	plate <- map[,c("Sample","Sample")]
	plate$region <- "v1v3"
	plate$study <- gsub("_.*","",map$Sample)
	write.table(plate,"plate",row.names=FALSE,quote=FALSE,col.names=FALSE,sep="\t");	
	'  | R --save
#load sumfiles
for sumfile in $(ls blast/*sum/*)
	do echo $sumfile;
	cat $sumfile   | awk '{printf $2"\t"$1"\t"$3"\t"$4"\n"}' > max_scores
	mysqlimport -u root -L $1 max_scores
	done


#load groups and plate
mysqlimport -u root -L $1 groups
mysqlimport -u root -L $1 plate


echo "UPDATE groups set region='v1v3' where region is null" | mysql -u root $1

#create max_scores_complete
echo "DROP TABLE IF EXISTS max_scores_complete;
	CREATE TABLE max_scores_complete select m.*,g.accession,g.barcode,g.region,p.sample,p.study,d.subset,d.otu,d.genus_group,d.phylum,d.class,d.genus from max_scores m left join groups g on m.query_id=g.accession left join plate p on g.barcode = p.barcode && g.region = p.region left join $blast_db d on m.db_id=d.accession; 
	ALTER TABLE max_scores_complete modify query_id VARCHAR(100) UNIQUE" | mysql -u root $1
