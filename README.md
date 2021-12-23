# Sanitize It Yourself
Source code for https://sanitizer.chemical.space/

## Usage
We recommend to use Docker.

1. Install Docker
2. Build the container image and run the image

```
git clone https://github.com/n-yoshikawa/molecule-sanitizer.git
cd molecule-sanitizer
docker build -t sanitizer .
docker run -p 8000:8000 -it sanitizer /bin/bash
```

3. Edit settings. You need to have keys and tokens for [Twitter API](https://developer.twitter.com/en/docs/authentication/oauth-1-0a) and [ORCID Public API](https://info.orcid.org/documentation/features/public-api/#easy-faq-2606) to run this application.  Callback URI should be set to `http://127.0.0.1:8000/complete/twitter/` if you run the server locally. If you don't need user authentication, consider using [local version](https://github.com/n-yoshikawa/molecule-sanitizer/tree/local) to avoid these settings.

```
cd molecule-sanitizer/
mv sanitizer/settings_local.example.py sanitizer/settings_local.py 
vim sanitizer/settings_local.py
```

4. Set up the database and run server
 
```
python3 manage.py makemigrations app
python3 manage.py migrate
python3 manage.py runserver 0.0.0.0:8000
```

You can access `http://127.0.0.1:8000/` in your browser.
