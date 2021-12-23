# Sanitize It Yourself (local version)
Source code for https://sanitizer.chemical.space/

## Usage
### Docker
We recommend to use Docker.

1. Install Docker
2. Build the container image and run the image

```
git clone https://github.com/n-yoshikawa/molecule-sanitizer.git
cd molecule-sanitizer
docker build -t sanitizer .
docker run -p 8000:8000 -it sanitizer /bin/bash
```

3. Edit settings. 

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

### Without Docker
If you don't want to use Docker, you can install dependencies by yourself.
```
pip install Django social-auth-app-django rdkit-pypi pandas annoy
```
