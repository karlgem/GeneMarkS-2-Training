
function getgroup() {
    local fnmod=$1;

    local group=$(grep TYPE $fnmod | cut -f2 -d "-");
    echo "$group";
}

function group_old_to_new () {
    local groupOld=$1;

    local groupNew;

    if [[ $groupOld == "a" ]]; then
        groupNew="d";
    if [[ $groupOld == "a2" ]]; then
        groupNew="d";
    elif [[ $groupOld == "b" ]]; then
        groupNew="c";
    elif [[ $groupOld == "c" ]]; then
        groupNew="b";
    elif [[ $groupOld == "d" ]]; then
        groupNew="a";
    elif [[ $groupOld == "e" ]]; then
        groupNew="x";
    fi
    
    echo "$groupNew";
}

function compare2gms2runs() {

    dir1=$1
    dir2=$2         # old version
    genomeList=$3

    cut -f1 $genomeList | while read -r genome; do

        path1=/storage3/w/gms2/data/$genome/runs/$dir1;
        path2=/storage3/w/gms2/data/$genome/runs/$dir2;

        # get group from 2 and convert to new naming convention
        local group2Old=$(getgroup $path2/GMS2.mod);
        local group2New=$(group_old_to_new $group2Old);

        local group1=$(getgroup $path1/GMS2.mod);

        echo -e "$group1\t$group2New\t$genome";
    done
}