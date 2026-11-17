#!/bin/bash

# collect static files
./manage.py collectstatic --noinput

# if DJANGO_DEBUG=1, run django dev server, else gunicorn
if [ "${DJANGO_DEBUG}" = "1" ]; then
	echo "Running Django development server..."
	./manage.py runserver 0.0.0:8000
else
	echo "Running Gunicorn server..."
	# run gunicorn server
	gunicorn \
		--worker-class asgi \
		--workers 4 --timeout 120 \
		--access-logfile - --error-logfile - \
		--bind 0.0.0.0:8000 \
		config.asgi:application
fi
