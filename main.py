import os

from flask import Flask, render_template, request, redirect, url_for, flash, session
from flask.json import jsonify
from flask_bootstrap import Bootstrap
from werkzeug.utils import secure_filename

from login import LoginForm
from register import RegisterForm
from session.session import Session
from user import User
from db import db
from werkzeug.security import generate_password_hash, check_password_hash
from flask_login import LoginManager, UserMixin, login_user, login_required, logout_user, current_user

from utils import profile
from utils.db_ops import checklogin
from utils.profile import new_dataset

app = Flask(__name__)
app.secret_key = 'amir'
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///username.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['UPLOAD_FOLDER'] = 'user_data'
bootstrap = Bootstrap(app)
login_manager = LoginManager()
login_manager.init_app(app)
login_manager.login_view = 'login'

APP_ROOT = os.path.dirname(os.path.abspath(__file__))

sess = Session(app)


@app.route('/')
def main():
    return redirect(url_for('login'))


@app.route('/upload', methods=['GET', 'POST'])
@login_required
def upload():
    sess.reset_user()
    user_dataset, user_configs, param_configs = profile.get_configs_files(APP_ROOT, session['user'])
    upload_page = '/upload_file_new_form.html' if user_configs else '/upload_file_first_form.html'

    if request.method == 'POST':
        if 'existingdataset' in request.form and request.form['existingdataset'] == 'on':
            pass # TODO
        else:
            new_dataset(request.form['datasetname'], request.files, APP_ROOT, session['user'], sess)
        return redirect(url_for('next'))

    if not user_configs:
        return render_template(upload_page, page=0)
    else:
        return render_template(upload_page, page=0, user_dataset=user_dataset, user_configs=user_configs)


@app.route('/dataset_configs')
def dataset_configs():
    user_dataset, user_configs, param_configs = profile.get_configs_files(APP_ROOT, session['user'])
    return jsonify(user_dataset=user_dataset, param_configs=param_configs, user_configs=user_configs)


@app.route('/next')
def next():
    return "<h1>hello world</h1>"


@login_manager.user_loader
def load_user(user_id):
    return User.query.get(int(user_id))


@app.route('/login', methods=['GET', 'POST'])
def login():
    form = LoginForm()
    if form.validate_on_submit():
        user = User.query.filter_by(username=form.username.data).first()
        if not checklogin(form, login_user, session, sess):
            return render_template('login.html', form=form, error='Invalid username or password')
        return redirect(url_for('upload'))
    return render_template('login.html', form=form)


@app.route('/signup', methods=['GET', 'POST'])
def signup():
    form = RegisterForm()
    if form.validate_on_submit():
        hashed_passwd = generate_password_hash(form.password.data, method='sha256')
        new_user = User(username=form.username.data, email=form.email.data, password=hashed_passwd)
        db.session.add(new_user)
        db.session.commit()
        return '<h1> user has been created!</h1>'
    return render_template('signup.html', form=form)


@app.route('/logout')
@login_required
def logout():
    logout_user()
    return redirect(url_for('index'))


def create_all():
    from db import db
    with app.app_context():
        db.init_app(app)
        db.create_all()
        db.session.commit()
        hashed_passwd = generate_password_hash('test12345', method='sha256')
        new_user = User(username='test', email='test@test.com', password=hashed_passwd)
        db.session.add(new_user)
        db.session.commit()
    return True


def create():
    from db import db
    with app.app_context():
        db.init_app(app)
        db.create_all()


db.init_app(app)

if __name__ == '__main__':
    app.run(debug=True)
