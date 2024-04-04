# Use an official Python runtime as a parent image
FROM python:3

# Set the working directory in the container
WORKDIR /tester


COPY ./requirements.txt /tester/requirements.txt



# Install any needed dependencies specified in requirements.txt
RUN pip install --no-cache-dir --upgrade -r requirements.txt

EXPOSE 80

COPY ./app /tester/app


CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "80"]


