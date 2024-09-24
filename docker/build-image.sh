
set -e

version=$1
harborHost=192.168.30.35:18082
repository=library
imageName=$harborHost/$repository/semd-builder:$version

if [ -z $1 ]; then
    echo "请传递版本号"
    exit 1  
fi

echo "Binding image....."
docker image build  -t $imageName -f docker/semd-builder-dockerfile docker

echo "Pushing image....."
docker push $imageName

echo "Congratulations! The semd-builder image:$imageName has been built successfully!!!"