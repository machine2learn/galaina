import os
import json
from forms.parameters_form import GeneralParameterForm

from flask import Flask, render_template, request, redirect, url_for, flash, session, send_file, send_from_directory
from flask.json import jsonify
from flask_bootstrap import Bootstrap
from werkzeug.utils import secure_filename
from threading import Thread
from login import LoginForm
from register import RegisterForm
from session.session import Session
from user import User
from db import db
from werkzeug.security import generate_password_hash, check_password_hash
from flask_login import LoginManager, UserMixin, login_user, login_required, logout_user, current_user
from flask_login import LoginManager, login_user, login_required, logout_user
from utils import profile
from utils.db_ops import checklogin
from utils.profile import new_dataset, set_dataset
from utils.parameters_util import check_causal_discovery_ob, set_form
from utils.sys_utils import delete_configs, delete_dataset
import python_r_pipeline

from config import config_reader  # FG import CustomMetaConfigParser


from werkzeug.datastructures import ImmutableMultiDict

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

path_to_meta_config = os.path.join(APP_ROOT, 'meta_config.ini')  # FG

sess = Session(app)


@app.route('/')
def main():
    return redirect(url_for('login'))


@app.route('/upload', methods=['GET', 'POST'])
@login_required
def upload():
    sess.reset_user()
    user_dataset, user_configs, param_configs = profile.get_configs_files(APP_ROOT, session['user'])
    if request.method == 'POST':
        if 'existingdataset' in request.form and request.form['existingdataset'] == 'on':
            set_dataset(APP_ROOT, session['user'], request.form['existing-select'],
                        json.loads(request.form['selected_config'])[0], sess)
        else:
            new_dataset(request.form['datasetname'], request.files, APP_ROOT, session['user'], sess)
        return redirect(url_for('parameters'))
    return render_template('/upload.html', page=0, user_dataset=user_dataset, user_configs=user_configs,
                           param_configs=param_configs)


@app.route('/dataset_configs')
def dataset_configs():
    user_dataset, user_configs, param_configs = profile.get_configs_files(APP_ROOT, session['user'])
    return jsonify(user_dataset=user_dataset, param_configs=param_configs, user_configs=user_configs)


@app.route('/parameters', methods=['GET', 'POST'])
def parameters():
    sess.open_log()
    form = GeneralParameterForm()
    if form.validate_on_submit():
        dict_parameters = request.form.to_dict()
        check_causal_discovery_ob(dict_parameters)
        sess.get_writer().load_all_sections()
        sess.get_writer().populate_config(dict_parameters)
        sess.get_writer().write_config()
        return redirect(url_for('run'))

    set_form(form, sess.get_writer())
    return render_template('parameters.html', form=form, page=1)


def threaded_function(writer, path):
    python_r_pipeline.run(writer, path)


@app.route('/run', methods=['GET', 'POST'])
def run():
    if request.method == 'POST':
        if 'remove_button' in request.form:
            _, dataset = sess.get_writer().get_info()
            delete_dataset(APP_ROOT, session['user'], dataset)
            return redirect(url_for('upload'))
        path, dataset_name = sess.get_writer().get_info()

        meta_config = config_reader.read_meta_config(path_to_meta_config)
        sess.get_writer().add_r_front_end(APP_ROOT, meta_config.get_path_r_binary_command())
        thread = Thread(
            target=threaded_function,
            args=(sess.get_writer().config, path)
        )
        thread.start()
        thread.join()
    return render_template('run.html', page=2)


@app.route("/download_pdf")
def download_pdf():
    pdf_path = sess.get_writer().get_output_path_fig()
    # FG better use os.path.basename
    filename = pdf_path.split('/')[-1]
    return send_file(pdf_path,
                     mimetype='text/csv',
                     attachment_filename=filename,
                     as_attachment=True)


@app.route("/show_pdf")
def show_pdf():
    pdf_path = sess.get_writer().get_output_path_fig()
    # FG better use os.path.basename
    filename = pdf_path.split('/')[-1]
    # FG better use os.path.dirname
    directory = '/'.join(pdf_path.split('/')[:-1])
    return send_from_directory(directory=directory,
                               mimetype='application/pdf',
                               filename=filename
                               )


@app.route('/stream')
@login_required
def stream():
    try:
        return jsonify(data=sess.read_log())
    except:
        return ''


@app.route('/delete_config', methods=['POST'])
@login_required
def delete_config():
    delete_configs(request.get_json()['config'], request.get_json()['dataset'], session['user'])
    user_dataset, user_configs, param_configs = profile.get_configs_files(APP_ROOT, session['user'])
    return jsonify(param_configs=param_configs, user_configs=user_configs)


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
def logout():
    logout_user()
    return redirect(url_for('login'))


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
    app.run(debug=True, host='127.0.0.1', port=5000)
