import os

from flask import Flask, render_template, request, redirect, url_for, flash
from flask_bootstrap import Bootstrap
from werkzeug.utils import secure_filename

from login import LoginForm
from register import RegisterForm
from user import User
from db import db
from werkzeug.security import generate_password_hash, check_password_hash
from flask_login import LoginManager, UserMixin, login_user, login_required, logout_user, current_user
from queue import Queue
from data import Data

from utils.profile import create_user_path

app = Flask(__name__)
app.secret_key = 'amir'
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///username.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['UPLOAD_FOLDER'] = 'user_data'
bootstrap = Bootstrap(app)
login_manager = LoginManager()
login_manager.init_app(app)
login_manager.login_view = 'login'
ALLOWED_EXTENSIONS = ['csv']

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/')
# @login_required
def index():
    return render_template('upload_file_new_form.html', page=0)




@app.route('/upload', methods=['POST'])
def upload():
    if request.method == 'POST':
        for filename, file in request.files.items():
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
        return redirect(url_for('next'))
    return redirect('/upload')


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
        if check_password_hash(user.password, form.password.data):
            login_user(user, remember=form.remember.data)
            create_user_path(user)
            return redirect(url_for('index'))
        return f'<h1>invalid username or password</h1>'
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


def create():
    from db import db
    with app.app_context():
        db.init_app(app)
        db.create_all()


if __name__ == '__main__':
    db.init_app(app)
    app.run(debug=True)
