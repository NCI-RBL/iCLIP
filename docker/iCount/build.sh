# Build image
docker build --no-cache -f Dockerfile --tag=ccbr_icount:v1.0.0 .

# Tag image with version and reset latest
docker tag ccbr_icount:v1.0.0 skchronicles/ccbr_icount:v1.0.0
docker tag ccbr_icount:v1.0.0 skchronicles/ccbr_icount
docker tag ccbr_icount:v1.0.0 nciccbr/ccbr_icount:v1.0.0
docker tag ccbr_icount:v1.0.0 nciccbr/ccbr_icount

# Push image to DockerHub
docker push skchronicles/ccbr_icount:v1.0.0
docker push skchronicles/ccbr_icount:latest
docker push nciccbr/ccbr_icount:v1.0.0
docker push nciccbr/ccbr_icount:latest
