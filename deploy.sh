#!/bin/bash
cd ~/project/

git config --global user.email "alcalan@fellows.iarc.fr"
git config --global user.name "Circle CI_$CIRCLE_PROJECT_REPONAME_$CIRCLE_BRANCH"
git add .
git status
git commit -m "Generated DAG [skip ci]"
git push origin $CIRCLE_BRANCH

curl -H "Content-Type: application/json" --data "{\"source_type\": \"Branch\", \"source_name\": \"$CIRCLE_BRANCH\"}" -X POST https://cloud.docker.com/api/build/v1/source/d3356607-1420-43f3-8955-8f9212639ba4/trigger/5f8ed292-0607-4279-88e7-e8d0c65df8e0/call/
