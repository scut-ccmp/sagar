PACKAGE="sagar"
VERSION_FILE=${PACKAGE}/version.py

version=$1
while true; do
  read -p "Release version ${version}? " yn
  case $yn in
    [Yy]* ) break;;
    [Nn]* ) exit;;
    * ) echo "Please answer yes or no.";;
  esac
done

# 字符串的1表示第一行
sed -i "1 s/__version__* =.*/__version__ = \"${version}\"/" $VERSION_FILE

# 得到当前所在分支
current_branch=`git rev-parse --abbrev-ref HEAD`
if [ "$current_branch" = "master" ]; then
  echo "ERROR: Do not modify and release at master branch."
  exit
fi

tag="v${version}"
relbranch="release-${version}"

echo Releasing version $version

# 创建并进入发布分支
git checkout -b $relbranch
git add ${VERSION_FILE}
git commit --no-verify -m "Release ${version}"

git tag -a $tag -m "Version $version"

# Merge into master

git checkout master
# --no-ff flag and git merge will always construct a merge instead of fast-forwarding.
# 结构上比较清晰
git merge --no-ff $relbranch

git checkout $current_branch
git merge --no-ff $relbranch

git branch -d $relbranch

# Push everything

# 需要current_branch 不是master 否则报错
# 要求开发者不会在master分支做修改，所有修改均来自其他分支。
git push --tags origin master $current_branch

# Release on pypi

rm -r dist
rm -r build
rm -r *.egg-info
python setup.py sdist
python setup.py bdist_wheel --universal

twine upload dist/*
